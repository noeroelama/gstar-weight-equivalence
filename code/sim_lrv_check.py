"""Diagnostic for Section 5.2: is the asymptotic variance right, and is the gap
between the estimated standard error and the Monte-Carlo standard deviation a
finite-sample bias of the long-run-variance (LRV) estimate?

Part 1. Population LRV. Simulate one very long series (T_LONG periods) from the
Part A DGP at the boundary of the equivalence region (|theta| = delta) for
phi11_mean, fit both candidate matrices, form the influence sequence psi_t, and
estimate its LRV with a Bartlett kernel and a wide bandwidth (and with the
Politis-White block bootstrap). The implied standard deviation at T=480 is
sqrt(LRV/480); compare it with sd_theta_hat in results/simA_size_power.csv.

Part 2. Large T. Replicate the T=2000 case R2 times and compare the mean MBB
standard error with the Monte-Carlo standard deviation.

Output: results/sim_lrv_check.txt
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gstarx as G          # noqa: E402
import inference as I       # noqa: E402
import dgp as D             # noqa: E402
from sim_study import pop_theta, lam_for_ratio   # noqa: E402

RES = HERE.parent / "results"
T_LONG = 400_000
T_BIG, R2, B2 = 2000, 500, 399
NAME = "phi11_mean"
SEED = 0


def main():
    lines = []
    log = lambda s: (print(s, flush=True), lines.append(s))
    dg = D.make_dgp(N=18, k=3, corr_range=0.6, seed=SEED)
    WA = dg["W0"]
    Wwrong = D.knn_from_coords(np.random.default_rng(SEED + 101).random((18, 2)), 3)
    th_far = abs(pop_theta(dg, WA, D.blend(WA, Wwrong, 1.0), NAME))
    delta = round(th_far / 2.2, 4)
    lam = lam_for_ratio(dg, WA, Wwrong, NAME, delta, 1.0)
    WB = D.blend(WA, Wwrong, lam)
    theta0 = pop_theta(dg, WA, WB, NAME)
    log(f"DGP as in Part A; functional={NAME}; delta={delta}; lam(ratio=1)={lam:.4f}; theta0={theta0:.4f}")

    # Part 1: population LRV from one long series
    t0 = time.time()
    rng = np.random.default_rng(123)
    Y = D.simulate(dg, T_LONG, rng)
    fA, fB = G.fit(Y, WA), G.fit(Y, WB)
    th, psi = G.theta_and_influence(fA, fB, NAME)
    n = psi.size
    for bw in [50, 100, 200, 400]:
        lrv = I.hac_lrv(psi, bw)
        log(f"[1] T_LONG={T_LONG}: theta_hat={th:.4f}; Bartlett LRV (bandwidth {bw}) = {lrv:.5f}; "
            f"implied sd at T=480 = {np.sqrt(lrv/480):.4f}; at T=240 = {np.sqrt(lrv/240):.4f}; "
            f"at T=120 = {np.sqrt(lrv/120):.4f}")
    Lpw = I.politis_white(psi)
    log(f"[1] Politis-White block length on the long series: {Lpw:.1f}  [{time.time()-t0:.0f}s]")

    simA = RES / "simA_size_power.csv"
    if simA.exists():
        df = pd.read_csv(simA)
        sub = df[(df.functional == NAME) & (np.isclose(df.ratio, 1.0))]
        for _, r in sub.iterrows():
            log(f"[1] Part A, ratio=1, T={int(r['T'])}: MC sd = {r['sd_theta_hat']:.4f}; "
                f"mean se MBB = {r['mean_se_mbb_pw']:.4f}; mean se HAC = {r['mean_se_hac']:.4f}")

    # Part 2: T = 2000
    t0 = time.time()
    ths, ses, ses_h = [], [], []
    for r in range(R2):
        rng = np.random.default_rng(SEED * 100003 + 50_000 + r)
        Yr = D.simulate(dg, T_BIG, rng)
        fAr, fBr = G.fit(Yr, WA), G.fit(Yr, WB)
        thr, psir = G.theta_and_influence(fAr, fBr, NAME)
        ths.append(thr)
        ses.append(I.mbb_se(psir, I.politis_white(psir), B2, rng))
        ses_h.append(I.hac_se(psir))
    ths, ses, ses_h = map(np.asarray, (ths, ses, ses_h))
    sd = ths.std(ddof=1)
    log(f"[2] T={T_BIG}, R={R2}, B={B2}: MC sd = {sd:.4f}; mean se MBB = {ses.mean():.4f} "
        f"({100*(ses.mean()/sd-1):+.1f}%); mean se HAC = {ses_h.mean():.4f} "
        f"({100*(ses_h.mean()/sd-1):+.1f}%)  [{time.time()-t0:.0f}s]")
    (RES / "sim_lrv_check.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
