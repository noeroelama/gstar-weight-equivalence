"""Per-replication records at the boundary |theta| = delta for Parts A, C and D.

The summary files of sim_study.py keep only means over replications. The mean of the
estimated standard error is dominated by a few replications whose fitted system is close
to the stability boundary (T = 120, long-run effects), so the summary tables of Section 5
report the median standard error next to a robust Monte-Carlo standard deviation. This
script recomputes the boundary cells with the same seeds as sim_study.py (so every
replication is identical) and saves one row per replication.
Output: ../results/sim_boundary_reps.csv
"""
from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import dgp as D                       # noqa: E402
import sim_study as S                 # noqa: E402


def designs(seed=0):
    dA = D.make_dgp(N=18, k=3, corr_range=0.6, seed=seed)
    wA = D.knn_from_coords(np.random.default_rng(seed + 101).random((18, 2)), 3)
    yield "A", dA, dA["W0"], wA, [120, 240, 480], seed * 100003
    dC = D.make_dgp(N=18, k=3, corr_range=0.6, seed=seed, graph="band")
    yield "C", dC, D.knn_from_coords(dC["coords"], 3), wA, [120, 240, 480], seed * 100003
    dD = D.make_dgp(N=40, k=5, corr_range=0.3, seed=seed + 40, graph="band", q=2, innov="t", df=5)
    wD = D.knn_from_coords(np.random.default_rng(seed + 140).random((40, 2)), 5)
    yield "D", dD, D.knn_from_coords(dD["coords"], 5), wD, [240, 480], seed * 100003


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", default="ACD")
    ap.add_argument("--R", type=int, default=1000)
    ap.add_argument("--B", type=int, default=399)
    a = ap.parse_args()
    out = S.RES / f"sim_boundary_reps_{a.parts}.csv"
    rows = []
    for part, dg, WA, Ww, Ts, base in designs():
        if part not in a.parts:
            continue
        for name in S.FUNCS:
            delta = round(abs(S.pop_theta(dg, WA, D.blend(WA, Ww, 1.0), name)) / 2.2, 4)
            lam = S.lam_for_ratio(dg, WA, Ww, name, delta, 1.0)
            WB = D.blend(WA, Ww, lam)
            theta0 = S.pop_theta(dg, WA, WB, name)
            for T in Ts:
                for r in range(a.R):
                    rec = S.one_rep(dg, WA, WB, name, T, theta0, delta, np.random.default_rng(base + r), a.B)
                    rec.update(part=part, functional=name, T=T, rep=r, theta0=theta0, delta=delta)
                    rows.append(rec)
                print(part, name, T, flush=True)
                pd.DataFrame(rows).to_csv(out, index=False)
