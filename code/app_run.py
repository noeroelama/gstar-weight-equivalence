"""Application: monthly rainfall over Indonesia (Section 6).

Steps
 1. Load panel (484 months x 40 locations), monthly anomalies.
 2. Build the four candidate W from the first m=240 months (no scored month used).
 3. In-sample GSTARX fits (anomalies, full sample, X lagged): coefficient summary,
    long-run and finite-horizon effects, spectral radius, residual diagnostics.
 4. Pairwise equivalence tests for phi11_mean, direct and indirect effects:
    theta_hat, se (vector MBB with Politis-White L; HAC), 90% interval, smallest
    delta at which equivalence is declared, decision at the pre-specified delta.
 5. Rolling-window one-step forecasts (m=240, M=244): RMSE per W (GSTAR and
    GSTARX), pairwise RMSE-difference equivalence tests, block-length sensitivity.
 6. Robustness: the same in-sample effects on raw (not seasonally adjusted) data,
    to show what the seasonal cycle does to the implied persistence.
Outputs in ../results/app_*.csv and ../results/app_summary.md
"""
from __future__ import annotations
import sys, time, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from data_io import load_panel, climatology, anomalies   # noqa: E402
import weights as Wm                                       # noqa: E402
import gstarx as G                                         # noqa: E402
import inference as I                                      # noqa: E402

RES = HERE.parent / "results"
RES.mkdir(exist_ok=True)
DATA = HERE.parent / "data"
M_WINDOW = 240
B = 1999
ALPHA = 0.05
REF = "geographic KNN"            # reference specification for relative margins
REL_MARGIN = 0.10                 # 10% of the reference value (Q1, Q2)
DELTA_RMSE = 1.0                  # mm per month (Q3) -- external tolerance, see paper
SEED = 2026


def ljung_box_reject(resid, lags=12, alpha=0.05):
    """Share of locations whose residuals reject white noise at `alpha`."""
    n, N = resid.shape
    rej = 0
    for i in range(N):
        e = resid[:, i] - resid[:, i].mean()
        acf = np.array([np.dot(e[:-k], e[k:]) / np.dot(e, e) for k in range(1, lags + 1)])
        Q = n * (n + 2) * np.sum(acf ** 2 / (n - np.arange(1, lags + 1)))
        rej += int(Q > chi2.ppf(1 - alpha, lags))
    return rej / N


def mean_abs_crosscorr(resid):
    C = np.corrcoef(resid.T)
    iu = np.triu_indices_from(C, 1)
    return float(np.abs(C[iu]).mean())


def main():
    t_all = time.time()
    rng = np.random.default_rng(SEED)
    P = load_panel()
    Y, X, month, idn = P["Y"], P["X"], P["month"], P["idn"]
    n, N = Y.shape
    log = []

    # ---- 2. candidate W from the first window
    rows0 = np.arange(M_WINDOW)
    Ytr = anomalies(Y[rows0], month[rows0], climatology(Y, month, rows0))
    Wset = Wm.build_all(Ytr, DATA, P["region_ids"])
    names = list(Wset)
    pd.DataFrame({k: [float((W > 0).sum(1).mean())] for k, W in Wset.items()},
                 index=["mean_nonzero_per_row"]).T.to_csv(RES / "app_W_summary.csv")

    # ---- 3. in-sample fits on full-sample anomalies
    clim_full = climatology(Y, month, np.arange(n))
    Ya = anomalies(Y, month, clim_full)
    fits = {k: G.fit(Ya, W, X) for k, W in Wset.items()}
    est_rows = []
    for k, f in fits.items():
        W = Wset[k]
        imp = G.impacts(f["beta"], W)
        row = dict(W=k, phi10_mean=f["beta"][:, 1].mean(), phi11_mean=G.phi11_mean(f["beta"]),
                   phi11_min=f["beta"][:, 2].min(), phi11_max=f["beta"][:, 2].max(),
                   n_phi11_negative=int((f["beta"][:, 2] < 0).sum()),
                   direct=imp["direct"], indirect=imp["indirect"], total=imp["total"], rho_A=imp["rho"],
                   ljung_box_reject_share=ljung_box_reject(f["resid"]),
                   mean_abs_resid_crosscorr=mean_abs_crosscorr(f["resid"]),
                   resid_sd_mm=float(f["resid"].std()))
        for h in (3, 6, 12, 24):
            row[f"indirect_h{h}"] = G.impacts(f["beta"], W, horizon=h)["indirect"]
        for name in ["phi11_mean", "direct", "indirect"]:
            ps = G.influence_single(f, G.FUNCTIONALS[name][1](f))
            row[f"se_{name}"] = I.mbb_se(ps, I.politis_white(ps), B, rng)
        est_rows.append(row)
    est = pd.DataFrame(est_rows)
    est.to_csv(RES / "app_insample_estimates.csv", index=False)
    log.append(f"in-sample fits done [{time.time()-t_all:.0f}s]")

    # ---- 4. pairwise equivalence for Q1/Q2 with a margin relative to the reference matrix:
    #      H1: |theta| < c tau(REF), tested as two one-sided tests on theta +/- c tau(REF)
    #      with their own influence sequences (inference.tost_relative). For the direct
    #      effect a second basis uses the W-dependent part, tau_dir(REF) - 1.
    ref = est.set_index("W")
    psi_ref = {name: G.influence_single(fits[REF], G.FUNCTIONALS[name][1](fits[REF]))
               for name in ["phi11_mean", "direct", "indirect"]}
    eq_rows = []
    z = I.norm.ppf(1 - ALPHA)
    for name in ["phi11_mean", "direct", "indirect"]:
        bases = [("ref", float(ref.loc[REF, name]))]
        if name == "direct":
            bases.append(("dyn", float(ref.loc[REF, name]) - 1.0))
        for a, b in itertools.combinations(names, 2):
            th, psi = G.theta_and_influence(fits[a], fits[b], name)
            row = dict(functional=name, W_A=a, W_B=b, theta_hat=th, se_hac=I.hac_se(psi))
            for tag, tau in bases:
                r = I.tost_relative(th, psi, tau, psi_ref[name], c=REL_MARGIN, B=B, rng=rng)
                rh = I.tost_relative(th, psi, tau, psi_ref[name], c=REL_MARGIN,
                                     V=I.hac_cov2(psi, psi_ref[name]))
                if tag == "ref":
                    row.update(se_mbb=r["se_theta"], L_pw=r["L_pw"], L_used=r["L_used"],
                               ci90_lo=th - z * r["se_theta"], ci90_hi=th + z * r["se_theta"],
                               z_se=z * r["se_theta"], min_delta_abs=abs(th) + z * r["se_theta"])
                row.update({f"tau_{tag}": tau, f"delta_{tag}": REL_MARGIN * tau,
                            f"c_min_{tag}": r["c_min_050"], f"equivalent_{tag}": r["equivalent_050"],
                            f"c_min_{tag}_a025": r["c_min_025"], f"equivalent_{tag}_a025": r["equivalent_025"],
                            f"c_min_{tag}_hac": rh["c_min_050"], f"equivalent_{tag}_hac": rh["equivalent_050"]})
            eq_rows.append(row)
    eq = pd.DataFrame(eq_rows)
    eq.to_csv(RES / "app_equivalence_Q12.csv", index=False)
    log.append(f"Q1/Q2 equivalence done [{time.time()-t_all:.0f}s]")

    # ---- 5. rolling forecasts
    fc = {}
    for use_x in (False, True):
        for k, W in Wset.items():
            fc[(k, use_x)] = G.rolling_forecast(Y, month, W, X if use_x else None, m=M_WINDOW, score_idx=idn)
    rmse_rows = [dict(W=k, model="GSTARX" if ux else "GSTAR", rmse_mm=v["rmse"], M=v["loss"].size)
                 for (k, ux), v in fc.items()]
    pd.DataFrame(rmse_rows).to_csv(RES / "app_rmse.csv", index=False)
    log.append(f"rolling forecasts done [{time.time()-t_all:.0f}s]")

    fc_rows = []
    for use_x in (False, True):
        for a, b in itertools.combinations(names, 2):
            la, lb = fc[(a, use_x)]["loss"], fc[(b, use_x)]["loss"]
            th, dt = I.rmse_diff_sequence(la, lb)
            L = I.politis_white(dt)
            se_mbb = I.mbb_se(dt, L, B, rng)
            se_hac = I.hac_se(dt)
            r = I.tost(th, se_mbb, DELTA_RMSE, ALPHA)
            r25 = I.tost(th, se_mbb, DELTA_RMSE, 0.025)
            row = dict(model="GSTARX" if use_x else "GSTAR", W_A=a, W_B=b, theta_rmse_mm=th,
                       delta_mm=DELTA_RMSE, L_pw=L, L_used=I.block_length(L), se_mbb=se_mbb, se_hac=se_hac,
                       equivalent_a025=r25["equivalent"], min_delta_mm_a025=r25["min_delta"],
                       ci90_lo=r["ci_lo"], ci90_hi=r["ci_hi"], equivalent_mbb=r["equivalent"],
                       equivalent_hac=I.tost(th, se_hac, DELTA_RMSE, ALPHA)["equivalent"],
                       p_tost_mbb=r["p_tost"], min_delta_mm=r["min_delta"])
            for Lfix in (2, 6, 12, 24):
                row[f"se_mbb_L{Lfix}"] = I.mbb_se(dt, Lfix, B, rng)
            fc_rows.append(row)
    pd.DataFrame(fc_rows).to_csv(RES / "app_equivalence_rmse.csv", index=False)
    log.append(f"RMSE equivalence done [{time.time()-t_all:.0f}s]")

    # ---- 6. robustness: raw data (no seasonal adjustment)
    raw_rows = []
    for k, W in Wset.items():
        f = G.fit(Y, W, X)
        imp = G.impacts(f["beta"], W)
        raw_rows.append(dict(W=k, phi10_mean=f["beta"][:, 1].mean(), phi11_mean=G.phi11_mean(f["beta"]),
                             direct=imp["direct"], indirect=imp["indirect"], rho_A=imp["rho"],
                             ljung_box_reject_share=ljung_box_reject(f["resid"])))
    pd.DataFrame(raw_rows).to_csv(RES / "app_insample_rawdata.csv", index=False)

    # ---- 7. data facts quoted in the paper: mean rainfall, seasonal pattern in the variance
    sd_c = np.array([[Ya[month == c, i].std(ddof=1) for i in range(N)] for c in range(1, 13)])   # (12, N)
    ratio = sd_c.max(axis=0) / sd_c.min(axis=0)
    facts = dict(mean_rain_mm=float(Y.mean()), mean_rain_scored_mm=float(Y[:, idn].mean()),
                 sd_ratio_median=float(np.median(ratio)), sd_ratio_min=float(ratio.min()),
                 sd_ratio_max=float(ratio.max()), sd_ratio_median_scored=float(np.median(ratio[idn])))
    pd.Series(facts).to_csv(RES / "app_data_facts.csv", header=["value"])

    # ---- summary
    with open(RES / "app_summary.md", "w", encoding="utf-8") as fh:
        fh.write(f"# Application results (generated {time.strftime('%Y-%m-%d %H:%M')})\n\n")
        fh.write(f"Panel {n} x {N}; scored {idn.size}; window m={M_WINDOW}; B={B}; alpha={ALPHA}; "
                 f"relative margin {REL_MARGIN} of {REF}; delta_RMSE={DELTA_RMSE} mm\n\n")
        fh.write("## In-sample estimates (anomalies, GSTARX)\n\n" + est.round(3).to_markdown(index=False) + "\n\n")
        fh.write("## Equivalence Q1/Q2\n\n" + eq.round(4).to_markdown(index=False) + "\n\n")
        fh.write("## RMSE\n\n" + pd.DataFrame(rmse_rows).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## Equivalence RMSE\n\n" + pd.DataFrame(fc_rows).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## Data facts\n\n" + pd.Series(facts).round(3).to_markdown() + "\n\n")
        fh.write("## Raw data (no seasonal adjustment)\n\n" + pd.DataFrame(raw_rows).round(3).to_markdown(index=False) + "\n\n")
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))
    print(f"total {time.time()-t_all:.0f}s")


if __name__ == "__main__":
    main()
