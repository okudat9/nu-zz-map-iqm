#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nu_zz_map_iqm.py  ——  IQM 共通 ν_ZZ 全マッピング（負対照内蔵・1ジョブ）
================================================================================
対応機体: Garnet(20q) / Emerald(54q) / 今後の IQM 機（coupling_map 自動取得）

設計の出自:
  ・符号付き差分Ramsey  ← nu_zz_tau_pinpoint_HW_20260613 / negctrl_SHORT_20260613
      <X>+i<Y> の複素位相差 Δφ=arg(z_on·conj(z_ref)) を τ で線形回帰。
      target の T1/T2・検出ずれは common-mode で相殺。符号が保存される。
  ・並列詰め込み        ← nu_zz_parallel_20260611
      互いに素なエッジを1回路に詰め、target ごとに別 clbit で読む。
  ・負対照内蔵          ← 今回の統合点（旧版は別スクリプトだった）
      「24/24 全部負」を二度と誤読しないため、非隣接ペアを同一ジョブに必ず入れる。

判定:
  非隣接 ≈ 0 かつ 隣接が有意   → 本物の直結 ZZ
  非隣接も隣接と同程度         → 差分prep共通系統（測定側の癖）。ここで停止
  両方 ≈ 0                     → この窓では解像せず（τ範囲かT2の問題）

実行（Windows）:
  1) 自己検証（クレジット消費なし・必ず先に実行）
       python nu_zz_map_iqm.py --sim
  2) 実機
       python nu_zz_map_iqm.py garnet
       python nu_zz_map_iqm.py emerald
     トークンは実行時に聞く（環境変数 IQM_TOKEN があればそれを使用）。

必要:
  pip install qiskit qiskit-iqm numpy scipy
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
# 機体プロファイル。T2 は既知較正の中央値（τ範囲の自動決定に使う）。
# 機体プロファイル。T2/T1 は IQM 公開較正の実測中央値（2026-07-20 取得）。
# 注: 出力のペア番号は 0 始まりのインデックス。IQM 表記では QB(index+1)。
#     例) 出力 (0, 1) = QB1-QB2
DEVICES = {
    "garnet":  dict(url="https://cocos.resonance.meetiqm.com/garnet",
                    n_qubits=20, t2_us=8.40,  t2_echo_us=19.48, t1_us=33.95),
    "emerald": dict(url="https://cocos.resonance.meetiqm.com/emerald",
                    n_qubits=54, t2_us=15.12, t2_echo_us=39.36, t1_us=52.03),
}

SHOTS        = 4096     # 削るなら 2048（コスト半分・誤差 1.3倍）
N_TAU        = 10       # τ点数。線形回帰なので密である必要はない
MAX_ADJACENT = 6        # 1ジョブに詰める隣接ペア数の上限
N_DISTANT    = 4        # 負対照ペア数（判定の基準線。3本以上を推奨）
TAU_SAFETY   = 2.2      # τ_max = TAU_SAFETY × T2。減衰込み sim で 2.0-2.5 が最良
RES_KHZ      = 0.6      # これ未満は「解像せず」扱い
MIN_HOPS     = 2        # 負対照の最小グラフ距離。2 = 直結カプラ無し（ZZは距離で指数減衰）

# --sim で注入する既知 ν_ZZ[kHz]。隣接だけに乗せ、非隣接は 0。
SIM_NU_ADJ   = [8.0, 3.0, 39.0, 12.0, 5.0, 20.0, 1.5, 30.0]
# ==============================================================================

SIMULATE = "--sim" in sys.argv
_pos = [a for a in sys.argv[1:] if not a.startswith("--")]
DEVICE = (_pos[0].lower() if _pos else "garnet")
if DEVICE not in DEVICES:
    raise SystemExit(f"未知の機体: {DEVICE}  （選択肢: {list(DEVICES)}）")
PROF = DEVICES[DEVICE]


# ------------------------------------------------------------------ backend --
def get_backend():
    """IQM Resonance に接続。--sim のときは None を返す。"""
    if SIMULATE:
        return None
    from iqm.qiskit_iqm import IQMProvider
    token = (os.environ.get("IQM_TOKEN", "").strip()
             or input("IQM Resonance token を貼って Enter: ").strip())
    if not token:
        raise SystemExit("token 未設定")
    provider = IQMProvider(PROF["url"], token=token)
    return provider.get_backend()


def get_coupling(backend):
    """無向エッジ集合を取得。--sim では線形鎖を仮定。"""
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
    """共有 qubit の無いエッジを貪欲に k 本選ぶ（1回路に同時に詰められる）。"""
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
    負対照: 直結していない (target, control) を選ぶ。
    adjacent で使った qubit は避け、グラフ距離が 3 hop 以上になる組を優先。
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
                continue                      # 直結はダメ
            if hops(t, c) < MIN_HOPS:
                continue                      # 近すぎるのもダメ
            if any(q in {t, c} for p in out for q in p):
                continue                      # 使い回さない
            out.append((t, c))
            if len(out) >= n_want:
                return out
    return out


# ------------------------------------------------------------------ circuit --
def build_parallel(pairs, tau_us, anc_state, basis, sim_nu=None):
    """
    全ペアを1回路に詰める。qr は物理 index に連番マップ。
      pair=(target, control)
      anc_state=1 → control を |1> に。delay 中に ZZ が target の位相を回す。
      basis 'X': H .. H        → P(0)=(1+cosφ)/2
      basis 'Y': H .. Sdg .. H → P(0)=(1+sinφ)/2
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
            # 物理ZZの模擬: ctrl=|1> のときだけ位相。減衰は測定後に別途適用。
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
    """SIM 用: T2 減衰で contrast を落とす（実機では自動的に起きる）。"""
    return exp_val * math.exp(-tau_us / t2_us)


def marginal_p0(counts, idx, n_cl):
    """clbit idx の P(=0)。qiskit の文字列は左が最上位 clbit。"""
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
    """加重最小二乗 y = a x + b。(a, b, var_a) を返す。"""
    x, y, w = map(np.asarray, (x, y, w))
    S, Sx, Sy = w.sum(), (w * x).sum(), (w * y).sum()
    Sxx, Sxy = (w * x * x).sum(), (w * x * y).sum()
    D = S * Sxx - Sx * Sx
    if abs(D) < 1e-300:
        return float("nan"), float("nan"), float("inf")
    return (S * Sxy - Sx * Sy) / D, (Sxx * Sy - Sx * Sxy) / D, S / D


def fit_nu(rows, shots):
    """
    1ペア分の点群から符号付き ν_ZZ[kHz] を抽出。
      z_ref = <X>+i<Y> (anc=0),  z_on = <X>+i<Y> (anc=1)
      Δφ(τ) = arg(z_on · conj(z_ref))  → unwrap → 傾き = 2π·ν_ZZ
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
    # 負対照を優先確保。足りなければ隣接を1本ずつ減らして再試行。
    n_adj = MAX_ADJACENT
    while n_adj >= 2:
        adjacent = greedy_disjoint(edges, n_adj)
        distant = pick_distant(edges, adjacent, N_DISTANT)
        if len(distant) >= N_DISTANT:
            break
        n_adj -= 1
    if len(distant) < 2:
        raise SystemExit("負対照を2本以上確保できません。MIN_HOPS を下げてください。")
    if n_adj < MAX_ADJACENT:
        print(f"  [調整] 負対照{N_DISTANT}本を確保するため隣接を "
              f"{MAX_ADJACENT}→{n_adj} に削減")
    pairs = adjacent + distant

    t2 = PROF["t2_us"]
    tau_max = TAU_SAFETY * t2
    taus = [round(x, 2) for x in np.linspace(2.0, tau_max, N_TAU)]
    sim_nu = ({p: SIM_NU_ADJ[i % len(SIM_NU_ADJ)] for i, p in enumerate(adjacent)}
              if SIMULATE else None)

    sched = schedule(pairs, taus)
    n_circ = len(sched)
    print("=" * 78)
    print(f"  ν_ZZ MAP — {'SIM' if SIMULATE else DEVICE.upper()}"
          f"  ({PROF['n_qubits']}q, T2={t2}µs)")
    print("=" * 78)
    print(f"  エッジ総数     : {len(edges)}")
    print(f"  隣接ペア       : {adjacent}")
    print(f"  非隣接（負対照）: {distant}")
    print(f"  τ              : {taus[0]} 〜 {taus[-1]} µs / {N_TAU}点"
          f"  (= {TAU_SAFETY}×T2)")
    print(f"  回路数         : {n_circ}  (1ジョブ)   shots={SHOTS}")
    pure = sum(s["tau_us"] for s in sched) * 1e-6 * SHOTS
    print(f"  純delay概算    : {pure:.2f}s  （+ゲート/読出/待機。実請求はこれより長い）")
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
    print(f"  {'pair':>10} {'(IQM)':>12} | {'種別':>6} | {'ν_ZZ (kHz)':>16} | {'contrast':>8}"
          + ("  | 真値(sim)" if SIMULATE else ""))
    print("  " + "-" * 74)
    rows = []
    for p in pairs:
        kind = "隣接" if p in adjacent else "非隣接"
        f = fits[p]
        val = f"{f['nu_khz']:+8.3f} ± {f['err_khz']:.3f}" if f else "n/a"
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
                         contrast=None if not f else round(f["contrast"], 4)))
    print("=" * 78)

    adj_raw = [fits[p]["nu_khz"] for p in adjacent if fits[p]]
    dis_raw = [fits[p]["nu_khz"] for p in distant if fits[p]]
    verdict = "INSUFFICIENT"
    if adj_raw and len(dis_raw) >= 2:
        # 負対照の「基準線」= 平均（オフセット）と広がり（ゆらぎ）
        base = float(np.mean(dis_raw))
        spread = float(np.std(dis_raw, ddof=1))
        cut = abs(base) + 2.0 * spread            # これを超えたら本物候補
        print(f"  負対照 基準線 : 平均 {base:+.3f} kHz   広がり ±{spread:.3f} kHz")
        print(f"  判定しきい値  : |ν_ZZ| > {cut:.3f} kHz  (= |平均| + 2×広がり)")
        print()
        real = [(p, fits[p]["nu_khz"]) for p in adjacent
                if fits[p] and abs(fits[p]["nu_khz"]) > cut]
        for p in adjacent:
            if not fits[p]:
                continue
            v = fits[p]["nu_khz"]
            mark = "★ 本物候補" if abs(v) > cut else "  基準線内（結合と判定できず）"
            print(f"    {str(p):>10} {v:+8.3f} kHz  {mark}")
        print()
        if abs(base) > 2.0 * spread and spread > 0:
            verdict = (f"【共通系統あり】負対照が0でなく {base:+.3f} kHz に偏っている。"
                       f"隣接の値からこのオフセットを差し引く必要がある。")
        elif real:
            verdict = (f"【本物の直結ZZ】{len(real)}/{len(adj_raw)} ペアが基準線を超えた。"
                       f"マップとして使える。")
        else:
            verdict = ("≈0 — 基準線を超える隣接ペアなし。"
                       "この機体/この窓では結合を解像できていない。")
        print(f"  → {verdict}")

    if SIMULATE:
        errs = [abs(fits[p]["nu_khz"] - sim_nu.get(p, 0.0))
                for p in pairs if fits[p]]
        typ_err = float(np.median([fits[p]["err_khz"] for p in pairs if fits[p]]))
        thr = max(3.0 * typ_err, 0.3)
        ok = max(errs) < thr if errs else False
        print(f"\n  self-test: {'PASS' if ok else 'FAIL'}"
              f"  (最大誤差 {max(errs):.3f} kHz / 許容 {thr:.3f} = 3σ)")
        print("  → PASS なら実機へ:  python nu_zz_map_iqm.py " + DEVICE)

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
