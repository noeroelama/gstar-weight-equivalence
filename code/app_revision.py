"""Regional effects and model-order diagnostics (Section 6.3 and 6.7 of the paper; Table S14).

Outputs in ../results/:
  rev_W.npz                 the four candidate W (from the first 240 months) and region ids
  rev_region_effects.csv    per-region phi11, direct effect M_ii, spillover received (row sum of M
                            off the diagonal) and emitted (column sum off the diagonal), per W
  rev_order_bic.csv         per W: share of locations whose BIC prefers GSTARX(1;1) over GSTARX(2;1),
                            and over the model without spatial lag (1;0); Ljung-Box rejection shares
  rev_order2_effects.csv    long-run average effects under GSTARX(2;1) (point estimates)
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from data_io import load_panel, climatology, anomalies   # noqa: E402
import weights as Wm                                       # noqa: E402
import gstarx as G                                         # noqa: E402
from app_run import ljung_box_reject, M_WINDOW             # noqa: E402

RES = HERE.parent / "results"
DATA = HERE.parent / "data"


def design_order(Y, W, X, p_own, p_sp, start):
    """Rows t = start..n-1; regressors: 1, Y_{t-l,i} (l<=p_own), (W Y_{t-l})_i (l<=p_sp), X_{t-1}."""
    n, N = Y.shape
    F = Y @ W.T
    rows = np.arange(start, n)
    cols = [np.ones((rows.size, N))]
    cols += [Y[rows - l] for l in range(1, p_own + 1)]
    cols += [F[rows - l] for l in range(1, p_sp + 1)]
    for j in range(X.shape[1]):
        cols.append(np.repeat(X[rows - 1, j][:, None], N, axis=1))
    return np.stack(cols, axis=2), Y[rows]


def bic(fit):
    n, k = fit["n"], fit["beta"].shape[1]
    rss = (fit["resid"] ** 2).sum(axis=0)
    return n * np.log(rss / n) + k * np.log(n)


def main():
    P = load_panel()
    Y, X, month = P["Y"], P["X"], P["month"]
    n, N = Y.shape
    rows0 = np.arange(M_WINDOW)
    Ytr = anomalies(Y[rows0], month[rows0], climatology(Y, month, rows0))
    Wset = Wm.build_all(Ytr, DATA, P["region_ids"])
    Ya = anomalies(Y, month, climatology(Y, month, np.arange(n)))
    np.savez(RES / "rev_W.npz", region_ids=np.array(P["region_ids"]),
             **{k.replace(" ", "_").replace("-", "_"): W for k, W in Wset.items()})

    reg, bic_rows, o2_rows = [], [], []
    for k, W in Wset.items():
        f = G.fit(Ya, W, X)
        imp = G.impacts(f["beta"], W)
        M = imp["M"]
        off = M - np.diag(np.diag(M))
        for i, rid in enumerate(P["region_ids"]):
            reg.append(dict(W=k, region_id=rid, country=P["countries"][i], phi11=f["beta"][i, 2],
                            direct=M[i, i], spill_received=off[i].sum(), spill_emitted=off[:, i].sum()))
        # order diagnostics on a common sample (t >= 2)
        f11 = G.ols(*design_order(Ya, W, X, 1, 1, 2))
        f21 = G.ols(*design_order(Ya, W, X, 2, 2, 2))
        f10 = G.ols(*design_order(Ya, W, X, 1, 0, 2))
        b11, b21, b10 = bic(f11), bic(f21), bic(f10)
        bic_rows.append(dict(W=k, share_bic_prefers_11_over_21=float((b11 <= b21).mean()),
                             share_bic_prefers_11_over_10=float((b11 <= b10).mean()),
                             lb_reject_11=ljung_box_reject(f11["resid"]),
                             lb_reject_21=ljung_box_reject(f21["resid"])))
        # GSTARX(2;1) long-run effects: M2 = (I - A1 - A2)^{-1}
        be = f21["beta"]
        A1 = np.diag(be[:, 1]) + np.diag(be[:, 3]) @ W
        A2 = np.diag(be[:, 2]) + np.diag(be[:, 4]) @ W
        Cmp = np.block([[A1, A2], [np.eye(N), np.zeros((N, N))]])
        M2 = np.linalg.inv(np.eye(N) - A1 - A2)
        tr = np.trace(M2)
        o2_rows.append(dict(W=k, direct=tr / N, indirect=(M2.sum() - tr) / N,
                            rho_companion=float(np.max(np.abs(np.linalg.eigvals(Cmp)))),
                            direct_order1=imp["direct"], indirect_order1=imp["indirect"]))
    pd.DataFrame(reg).to_csv(RES / "rev_region_effects.csv", index=False)
    pd.DataFrame(bic_rows).to_csv(RES / "rev_order_bic.csv", index=False)
    pd.DataFrame(o2_rows).to_csv(RES / "rev_order2_effects.csv", index=False)
    print(pd.DataFrame(bic_rows).round(3).to_string())
    print(pd.DataFrame(o2_rows).round(3).to_string())


if __name__ == "__main__":
    main()
