#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nu_zz_map_iqm.py  --  Residual ZZ mapping for IQM processors (one job, with a
                      built-in negative control)
================================================================================
Devices: Garnet (20q), Emerald (54q), and any future IQM backend -- the coupling
map, pair selection and delay range are all derived from the backend.

Design:
  * Signed differential Ramsey
      Linear fit of the complex phase difference dphi = arg(z_on * conj(z_ref))
      against tau, where z = <X> + i<Y>. Multiplicative corruption (target
      dephasing, contrast loss) cancels in the ratio. The sign is preserved.
  * Parallel packing
      Disjoint edges share one circuit; each target is read on its own clbit.
  * Built-in negative control
      Non-adjacent pairs are measured in the same job, every time. They set the
      baseline against which the real edges are judged.

Verdicts:
  controls near zero, adjacent significant  -> real direct ZZ coupling
  controls as large as adjacent             -> common systematic; stop here
  both near zero                            -> unresolved in this window

Usage (Windows):
  1) Self-test first. Costs nothing.
       python nu_zz_map_iqm.py --sim garnet
       python nu_zz_map_iqm.py --sim emerald
  2) Hardware
       python nu_zz_map_iqm.py garnet
       python nu_zz_map_iqm.py emerald
     The token is requested at runtime, or read from IQM_TOKEN.

Requires:
  pip install qiskit qiskit-iqm numpy scipy

Author : Takeshi Okuda  (ORCID 0009-0006-7449-202X)
License: Apache 2.0 -- see LICENSE
Repo   : https://github.com/okudat9/nu-zz-map-iqm
================================================================================
"""

import os
import sys
import json
import math
import datetime as _dt
from pathlib import Path

import numpy as np

# ============================ CONFIG ==========================================
# Device profiles (median calibration values, used to set the tau range).
# Device profiles. T2/T1 are median values taken from IQM's published
# calibration data (retrieved 2026-07-20).
# NOTE: pairs are reported as zero-based coupling-map indices.
#       IQM labels qubits from 1, so output (0, 1) means QB1-QB2.
DEVICES = {
    "garnet":  dict(url="https://cocos.resonance.meetiqm.com/garnet",
                    n_qubits=20, t2_us=8.40,  t2_echo_us=19.48, t1_us=33.95),
    "emerald": dict(url="https://cocos.resonance.meetiqm.com/emerald",
                    n_qubits=54, t2_us=15.12, t2_echo_us=39.36, t1_us=52.03),
}

SHOTS        = 4096     # halve to 2048 for half the cost, ~1.3x larger error
N_TAU        = 10       # number of tau points; a linear fit does not need many
MAX_ADJACENT = 6        # max coupled edges packed into one job
N_DISTANT    = 4        # negative controls; they set the baseline. 3+ recommended
TAU_SAFETY   = 2.2      # tau_max = TAU_SAFETY x T2. With decay included the
                        # optimum is flat from 2.0 to 2.5; error rises outside
RES_KHZ      = 0.6      # below this magnitude, treated as unresolved
MIN_HOPS     = 2        # min graph distance for controls. 2 = no direct coupler
                        # (ZZ falls off exponentially with distance)

# Known nu_ZZ [kHz] injected in --sim mode, applied to adjacent pairs only.
# Controls get zero, so the self-test checks both signal and baseline.
SIM_NU_ADJ   = [8.0, 3.0, 39.0, 12.0, 5.0, 20.0, 1.5, 30.0]
# ==============================================================================

SIMULATE = "--sim" in sys.argv
_pos = [a for a in sys.argv[1:] if not a.startswith("--")]
DEVICE = (_pos[0].lower() if _pos else "garnet")
if DEVICE not in DEVICES:
    raise SystemExit(f"unknown device: {DEVICE}  (choose from {list(DEVICES)})")
PROF = DEVICES[DEVICE]


# ------------------------------------------------------------------ backend --
def get_backend():
    """Connect to IQM Resonance. Returns None in --sim mode."""
    if SIMULATE:
        return None
    from iqm.qiskit_iqm import IQMProvider
    token = (os.environ.get("IQM_TOKEN", "").strip()
             or input("Paste your IQM Resonance token and press Enter: ").strip())
    if not token:
        raise SystemExit("no token provided")
    provider = IQMProvider(PROF["url"], token=token)
    return provider.get_backend()


def get_coupling(backend):
    """Fetch the undirected edge set. In --sim mode a linear chain is assumed."""
    if SIMULATE:
        n = PROF["n_qubits"]
        return sorted({(i, i + 1) for i in range(n - 1)})
    try:
        edges = list(backend.coupling_map.get_edges())
    except Exception:
        edges = backend.configuration().coupling_map
    return sorted({tuple(sorted(e)) for e in edges})


# ------------------------------------------------------------ pair selection --
def greedy_disjoint(edges, k):
    """Greedily pick k edges sharing no qubits, so they can run in parallel."""
    used, out = set(), []
    for a, b in edges:
        if a in used or b in used:
            continue
        out.append((a, b)); used.update((a, b))
        if len(out) >= k:
            break
    return out


def pick_distant(edges, adjacent, n_want):
    """
    Negative controls: pick (target, control) pairs with no direct coupler.
    Qubits already used by adjacent pairs are excluded; candidates must be at
    least MIN_HOPS apart on the coupling graph.
    """
    nodes = sorted({q for e in edges for q in e})
    adj = {q: set() for q in nodes}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)

    def hops(s, t, cap=4):
        seen, frontier = {s}, [s]
        for d in range(1, cap + 1):
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v == t:
                        return d
                    if v not in seen:
                        seen.add(v); nxt.append(v)
            frontier = nxt
            if not frontier:
                break
        return cap + 1

    busy = {q for p in adjacent for q in p}
    free = [q for q in nodes if q not in busy]
    out = []
    for i, t in enumerate(free):
        for c in free[i + 1:]:
            if c in adj[t]:
                continue                      # directly coupled -- not a control
            if hops(t, c) < MIN_HOPS:
                continue                      # too close on the graph
            if any(q in {t, c} for p in out for q in p):
                continue                      # qubit already used
            out.append((t, c))
            if len(out) >= n_want:
                return out
    return out


# ------------------------------------------------------------------ circuit --
def build_parallel(pairs, tau_us, anc_state, basis, sim_nu=None):
    """
    Pack every pair into one circuit; physical qubits are remapped to a
    contiguous register.
      pair = (target, control)
      anc_state=1 -> control prepared in |1>, so ZZ rotates the target phase
                     during the delay
      basis 'X': H .. H        -> P(0) = (1 + cos phi) / 2
      basis 'Y': H .. Sdg .. H -> P(0) = (1 + sin phi) / 2
    """
    from qiskit import QuantumCircuit
    qubits = sorted({q for p in pairs for q in p})
    qmap = {q: i for i, q in enumerate(qubits)}
    qc = QuantumCircuit(len(qubits), len(pairs))
    dt_ns = round(tau_us * 1000)

    for i, (tgt, ctrl) in enumerate(pairs):
        t, c = qmap[tgt], qmap[ctrl]
        if anc_state == 1:
            qc.x(c)
        qc.h(t)
        if SIMULATE:
            # Simulate physical ZZ: phase only when ctrl=|1>. Decay applied later.
            nu = (sim_nu or {}).get((tgt, ctrl), 0.0)
            if anc_state == 1 and nu:
                qc.rz(2 * math.pi * nu * 1e3 * tau_us * 1e-6, t)
        else:
            qc.delay(dt_ns, t, unit="ns")
            qc.delay(dt_ns, c, unit="ns")
        if basis == "Y":
            qc.sdg(t)
        qc.h(t)
        qc.measure(t, i)
    return qc


def schedule(pairs, taus):
    return [dict(tau_us=t, anc=a, basis=b)
            for a in (0, 1) for b in ("X", "Y") for t in taus]


# -------------------------------------------------------------------- runner --
def run(circuits, backend):
    if SIMULATE:
        from qiskit_aer import AerSimulator
        from qiskit import transpile
        sim = AerSimulator()
        res = sim.run(transpile(circuits, sim), shots=SHOTS).result()
        return [res.get_counts(i) for i in range(len(circuits))], "SIM"
    from iqm.qiskit_iqm import transpile_to_IQM
    isa = [transpile_to_IQM(c, backend) for c in circuits]
    job = backend.run(isa, shots=SHOTS)
    jid = getattr(job, "job_id", lambda: "n/a")()
    print(f"  job_id = {jid}")
    res = job.result()
    return [res.get_counts(i) for i in range(len(circuits))], jid


def apply_decay(exp_val, tau_us, t2_us):
    """Apply T2 decay in --sim mode; on hardware this happens by itself."""
    return exp_val * math.exp(-tau_us / t2_us)


def marginal_p0(counts, idx, n_cl):
    """P(0) for classical bit idx. In Qiskit bitstrings the leftmost
    character is the highest clbit."""
    tot = sum(counts.values())
    if not tot:
        return 0.5
    z = 0
    for b, c in counts.items():
        s = b.replace(" ", "")
        if s[n_cl - 1 - idx] == "0":
            z += c
    return z / tot


# ---------------------------------------------------------------- analysis ---
def wls(x, y, w):
    """Weighted least squares for y = a x + b. Returns (a, b, var_a)."""
    x, y, w = map(np.asarray, (x, y, w))
    S, Sx, Sy = w.sum(), (w * x).sum(), (w * y).sum()
    Sxx, Sxy = (w * x * x).sum(), (w * x * y).sum()
    D = S * Sxx - Sx * Sx
    if abs(D) < 1e-300:
        return float("nan"), float("nan"), float("inf")
    return (S * Sxy - Sx * Sy) / D, (Sxx * Sy - Sx * Sxy) / D, S / D


def fit_nu(rows, shots):
    """
    Extract signed nu_ZZ [kHz] for one pair.
      z_ref = <X> + i<Y>   (control in |0>)
      z_on  = <X> + i<Y>   (control in |1>)
      dphi(tau) = arg(z_on * conj(z_ref)) -> unwrap -> slope = 2*pi*nu_ZZ
    Multiplicative corruption (target dephasing, contrast loss) cancels in the
    ratio; additive corruption (asymmetric readout) does not. See README.
    """
    by_tau = {}
    for r in rows:
        by_tau.setdefault(r["tau_us"], {})[(r["anc"], r["basis"])] = r["exp"]

    ts, dphi, wts, contrast = [], [], [], []
    for t in sorted(by_tau):
        d = by_tau[t]
        if len(d) < 4:
            continue
        Xr, Yr = d[(0, "X")], d[(0, "Y")]
        Xo, Yo = d[(1, "X")], d[(1, "Y")]
        zr, zo = complex(Xr, Yr), complex(Xo, Yo)
        rr, ro = abs(zr), abs(zo)
        if rr < 1e-6 or ro < 1e-6:
            continue

        def sig(X, Y, r):
            return math.sqrt(max(1 - X * X, 1e-6) / shots
                             + max(1 - Y * Y, 1e-6) / shots) / r
        s = math.hypot(sig(Xr, Yr, rr), sig(Xo, Yo, ro))
        ts.append(t * 1e-6)
        dphi.append(np.angle(zo * np.conj(zr)))
        wts.append(1.0 / max(s, 1e-9) ** 2)
        contrast.append(ro)

    if len(ts) < 3:
        return None
    a, _, var = wls(ts, np.unwrap(np.array(dphi)), wts)
    return dict(nu_khz=a / (2 * np.pi) / 1e3,
                err_khz=math.sqrt(var) / (2 * np.pi) / 1e3,
                n=len(ts), contrast=float(np.mean(contrast)))


# ------------------------------------------------------------------- main ----
def main():
    backend = get_backend()
    edges = get_coupling(backend)
    # Secure the controls first; drop one adjacent edge at a time if needed.
    n_adj = MAX_ADJACENT
    while n_adj >= 2:
        adjacent = greedy_disjoint(edges, n_adj)
        distant = pick_distant(edges, adjacent, N_DISTANT)
        if len(distant) >= N_DISTANT:
            break
        n_adj -= 1
    if len(distant) < 2:
        raise SystemExit("cannot secure 2+ negative controls. Lower MIN_HOPS.")
    if n_adj < MAX_ADJACENT:
        print(f"  [adjusted] adjacent edges reduced {MAX_ADJACENT} -> {n_adj} "
              f"to secure {N_DISTANT} controls")
    pairs = adjacent + distant

    t2 = PROF["t2_us"]
    tau_max = TAU_SAFETY * t2
    taus = [round(x, 2) for x in np.linspace(2.0, tau_max, N_TAU)]
    sim_nu = ({p: SIM_NU_ADJ[i % len(SIM_NU_ADJ)] for i, p in enumerate(adjacent)}
              if SIMULATE else None)

    sched = schedule(pairs, taus)
    n_circ = len(sched)
    print("=" * 78)
    print(f"  RESIDUAL ZZ MAP — {'SIM' if SIMULATE else DEVICE.upper()}"
          f"  ({PROF['n_qubits']}q, T2={t2} us)")
    print("=" * 78)
    print(f"  edges on chip    : {len(edges)}")
    print(f"  adjacent pairs   : {adjacent}")
    print(f"  controls         : {distant}")
    print(f"  tau              : {taus[0]} - {taus[-1]} us / {N_TAU} points"
          f"  (= {TAU_SAFETY} x T2)")
    print(f"  circuits         : {n_circ}  (one job)   shots={SHOTS}")
    pure = sum(s["tau_us"] for s in sched) * 1e-6 * SHOTS
    print(f"  pure delay       : {pure:.2f}s  (plus gates/readout/queue; the\n                     billed time will be longer)")
    print()

    circs = [build_parallel(pairs, s["tau_us"], s["anc"], s["basis"], sim_nu)
             for s in sched]
    counts, jid = run(circs, backend)

    n_cl = len(pairs)
    for s, c in zip(sched, counts):
        s["_counts"] = c
    per_pair = {p: [] for p in pairs}
    for i, p in enumerate(pairs):
        for s in sched:
            e = 2.0 * marginal_p0(s["_counts"], i, n_cl) - 1.0
            if SIMULATE:
                e = apply_decay(e, s["tau_us"], t2)
            per_pair[p].append(dict(tau_us=s["tau_us"], anc=s["anc"],
                                    basis=s["basis"], exp=e))

    fits = {p: fit_nu(per_pair[p], SHOTS) for p in pairs}

    print("=" * 78)
    print(f"  {'pair':>10} {'(IQM)':>12} | {'type':>8} | {'nu_ZZ (kHz)':>16} | {'contrast':>8}"
          + ("  | true (sim)" if SIMULATE else ""))
    print("  " + "-" * 74)
    rows = []
    for p in pairs:
        kind = "adjacent" if p in adjacent else "control"
        f = fits[p]
        val = f"{f['nu_khz']:+8.3f} +/- {f['err_khz']:.3f}" if f else "n/a"
        ct = f"{f['contrast']:.3f}" if f else "n/a"
        extra = ""
        if SIMULATE:
            extra = f"  | {sim_nu.get(p, 0.0):.1f}"
        qb = f"QB{p[0]+1}-QB{p[1]+1}"
        print(f"  {str(p):>10} {qb:>12} | {kind:>6} | {val:>16} | {ct:>8}{extra}")
        rows.append(dict(pair=list(p), iqm_label=qb,
                         kind=("adjacent" if p in adjacent else "distant"),
                         nu_khz=None if not f else round(f["nu_khz"], 4),
                         err_khz=None if not f else round(f["err_khz"], 4),
                         contrast=None if not f else round(f["contrast"], 4),
                         sigma=None if not (f and "sigma" in f) else round(f["sigma"], 2),
                         resolved=None if not (f and "sigma" in f) else bool(f["sigma"] > 2.0)))
    print("=" * 78)

    adj_raw = [fits[p]["nu_khz"] for p in adjacent if fits[p]]
    dis_raw = [fits[p]["nu_khz"] for p in distant if fits[p]]
    verdict = "INSUFFICIENT"
    if adj_raw and len(dis_raw) >= 2:
        # Control baseline, tested against the MEASUREMENT errors rather than
        # the scatter of a few points -- a handful of values can cluster by
        # chance and make a pure-noise mean look significant.
        dv = np.array(dis_raw)
        de = np.array([fits[p]["err_khz"] for p in distant if fits[p]])
        wv = 1.0 / de**2
        base = float((dv * wv).sum() / wv.sum())          # weighted mean
        base_se = float(1.0 / np.sqrt(wv.sum()))          # its standard error
        sig = abs(base) / base_se if base_se > 0 else 0.0
        print(f"  control baseline : {base:+.3f} +/- {base_se:.3f} kHz "
              f"({sig:.2f} sigma from zero)")
        print(f"  test             : each pair vs baseline, using its own error")
        print()

        # Per-pair significance. A pair with a large error bar can sit far from
        # the baseline and still be unresolved -- one global threshold hides
        # that, so each pair is judged against its own precision.
        real = []
        for p in adjacent:
            f = fits[p]
            if not f:
                continue
            v = f["nu_khz"]
            dev = abs(v - base)
            unc = math.hypot(f["err_khz"], base_se)
            z = dev / unc if unc > 0 else 0.0
            f["sigma"] = z
            ok = z > 2.0
            if ok:
                real.append((p, v))
            mark = (f"* coupled          ({z:.1f} sigma)" if ok
                    else f"  not distinguishable ({z:.1f} sigma)")
            print(f"    {str(p):>10} {v:+8.3f} +/- {f['err_khz']:.3f} kHz  {mark}")
        print()
        if sig > 2.0:
            verdict = (f"COMMON SYSTEMATIC: controls sit {base:+.3f} kHz from zero "
                       f"({sig:.1f} sigma). Subtract this offset from the adjacent "
                       f"values, or fix the state preparation.")
        elif real:
            verdict = (f"REAL DIRECT ZZ: {len(real)}/{len(adj_raw)} adjacent pairs "
                       f"are resolved above the control baseline. "
                       f"{'All' if len(real) == len(adj_raw) else 'The remainder'} "
                       f"{'are usable.' if len(real) == len(adj_raw) else 'need more shots or a longer T2.'}")
        else:
            verdict = ("UNRESOLVED: no adjacent pair is distinguishable from the "
                       "control baseline. Raise the shots, or the coupling is "
                       "below the resolution of this device.")
        print(f"  -> {verdict}")

    if SIMULATE:
        errs = [abs(fits[p]["nu_khz"] - sim_nu.get(p, 0.0))
                for p in pairs if fits[p]]
        typ_err = float(np.median([fits[p]["err_khz"] for p in pairs if fits[p]]))
        thr = max(3.0 * typ_err, 0.3)
        ok = max(errs) < thr if errs else False
        print(f"\n  self-test: {'PASS' if ok else 'FAIL'}"
              f"  (max error {max(errs):.3f} kHz / tolerance {thr:.3f} = 3 sigma)")
        print("  -> if PASS, run on hardware:  python nu_zz_map_iqm.py " + DEVICE)

    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = Path(__file__).resolve().parent / "results" / ts[:8]
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"nu_zz_map_{'SIM' if SIMULATE else DEVICE}_{ts}.json"
    out.write_text(json.dumps(dict(
        tool="nu_zz_map_iqm", device=("SIM" if SIMULATE else DEVICE),
        job_id=jid, shots=SHOTS, taus_us=taus, t2_us=t2,
        n_edges=len(edges),
        adjacent=[list(p) for p in adjacent],
        distant=[list(p) for p in distant],
        results=rows, verdict=verdict, generated_at=ts),
        ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n  saved: {out}")


if __name__ == "__main__":
    main()
