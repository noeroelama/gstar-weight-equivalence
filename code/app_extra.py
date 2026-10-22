"""Additional application analyses (robustness checks in Section 6 and the Supplementary Material).

 A. Uniform (equal-weight) W as a fifth candidate; pairwise equivalence vs the four.
 B. Common-factor check: add the lagged mean anomaly of the other N-1 regions,
    (W_U Y_{t-1})_i with W_U the uniform matrix (zero diagonal), as a regressor; recompute
    effects with A = diag(phi10) + diag(phi11) W + diag(phi_m) W_U; pairwise equivalence
    of direct/indirect effects with margins relative to the geographic-KNN value.
 C. Similarity of the spatial-lag series W Y across candidates (correlations).
 D. J-test (Davidson-MacKinnon form): per-location augmented regression with HAC t,
    share of locations rejecting; pooled rho with vector-MBB standard error.
 E. Benchmark RMSEs: window climatology only; per-location AR(1) + climate (no W).
 F. Sensitivity: KNN k in {3, 5, 8} (from centroids, inverse-distance weights within the
   k neighbours, so k=5 reproduces the main geographic KNN matrix); window m in {180, 240, 300}.
 G. In-sample fits on months 241-484 only (disjoint from the W-construction sample).
 H. Dynamic-part margin for the direct effect (direct - 1 under KNN).
 I. Standardised anomalies (divided by the calendar-month standard deviation), to check
    that the seasonal pattern in the variance does not drive the results.
All relative margins use inference.tost_relative (two one-sided tests on theta +/- c tau_ref).
Outputs: ../results/extra_*.csv and ../results/extra_summary.md
"""
from __future__ import annotations
import sys, time, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from data_io import load_panel, climatology, anomalies   # noqa: E402
import weights as Wm                                       # noqa: E402
import gstarx as G                                         # noqa: E402
import inference as I                                      # noqa: E402

RES = HERE.parent / "results"; DATA = HERE.parent / "data"
M_WINDOW, B, ALPHA, SEED = 240, 1999, 0.05, 2026
Z1 = norm.ppf(1 - ALPHA)
REF = "geographic KNN"; REL = 0.10

# ------------------------------------------------------------ generalised effects


def transition_multi(beta, Ws, cols):
    """A = diag(beta[:,1]) + sum_k diag(beta[:, cols[k]]) Ws[k]."""
    A = np.diag(beta[:, 1])
    for W, c in zip(Ws, cols):
        A = A + np.diag(beta[:, c]) @ W
    return A


def impacts_multi(beta, Ws, cols):
    A = transition_multi(beta, Ws, cols); N = A.shape[0]
    M = np.linalg.solve(np.eye(N) - A, np.eye(N))
    tr, tot = np.trace(M), M.sum()
    return dict(direct=tr / N, indirect=(tot - tr) / N, rho=float(np.max(np.abs(np.linalg.eigvals(A)))), M=M)


def impact_grad_multi(beta, Ws, cols, which):
    A = transition_multi(beta, Ws, cols); N = A.shape[0]
    M = np.linalg.solve(np.eye(N) - A, np.eye(N)); M2 = M @ M
    rowsum, colsum = M.sum(axis=0), M.sum(axis=1)
    Gr = np.zeros_like(beta)
    dtr = {1: np.diag(M2)}; dtot = {1: rowsum * colsum}
    for W, c in zip(Ws, cols):
        dtr[c] = np.diag(W @ M2); dtot[c] = rowsum * (W @ colsum)
    for c in dtr:
        if which == "direct":
            Gr[:, c] = dtr[c] / N
        elif which == "indirect":
            Gr[:, c] = (dtot[c] - dtr[c]) / N
    return Gr


def fit_general(Ya, regs, X=None):
    """OLS with regressors: 1, Y_{t-1,i}, then each (R Y_{t-1})_i for R in regs, then X_{t-1}."""
    n, N = Ya.shape
    cols = [np.ones((n - 1, N)), Ya[:-1]]
    for R in regs:
        cols.append((Ya @ R.T)[:-1])
    if X is not None:
        for j in range(X.shape[1]):
            cols.append(np.repeat(X[:-1, j][:, None], N, axis=1))
    Z = np.stack(cols, axis=2)
    out = G.ols(Z, Ya[1:]); out["Z"] = Z
    return out


def theta_psi(fA, fB, WsA, WsB, cols, which):
    hA = impacts_multi(fA["beta"], WsA, cols)[which]; hB = impacts_multi(fB["beta"], WsB, cols)[which]
    gA = impact_grad_multi(fA["beta"], WsA, cols, which); gB = impact_grad_multi(fB["beta"], WsB, cols, which)
    return hA - hB, G.influence_single(fA, gA) - G.influence_single(fB, gB)


def rel_row(th, psi, tau_ref, psi_ref, rng, prefix=""):
    """Relative-margin equivalence (|theta| < c tau_ref) with the vector MBB, plus the HAC
    decision and the plain 90% interval for theta."""
    r = I.tost_relative(th, psi, tau_ref, psi_ref, c=REL, B=B, rng=rng)
    rh = I.tost_relative(th, psi, tau_ref, psi_ref, c=REL, V=I.hac_cov2(psi, psi_ref))
    out = dict(theta_hat=th, se_mbb=r["se_theta"], L_pw=r["L_pw"], ci90_lo=th - Z1 * r["se_theta"],
               ci90_hi=th + Z1 * r["se_theta"], z_se=Z1 * r["se_theta"])
    out.update({f"{prefix}tau_ref": tau_ref, f"{prefix}c_min": r["c_min_050"],
                f"{prefix}equivalent": r["equivalent_050"], f"{prefix}c_min_a025": r["c_min_025"],
                f"{prefix}equivalent_a025": r["equivalent_025"], f"{prefix}equivalent_hac": rh["equivalent_050"]})
    return out


def infl(f, W, which):
    return G.influence_single(f, G.impact_grad(f["beta"], W, which))


def eq_row(th, psi, delta, rng):
    L = I.politis_white(psi); se = I.mbb_se(psi, L, B, rng); r = I.tost(th, se, delta, ALPHA)
    return dict(theta_hat=th, se_mbb=se, se_hac=I.hac_se(psi), L_pw=L, ci90_lo=r["ci_lo"], ci90_hi=r["ci_hi"],
                delta=delta, equivalent=r["equivalent"], min_delta=r["min_delta"], z_se=Z1 * se)


def main():
    t0 = time.time(); rng = np.random.default_rng(SEED)
    P = load_panel(); Y, X, month, idn = P["Y"], P["X"], P["month"], P["idn"]; n, N = Y.shape
    rows0 = np.arange(M_WINDOW)
    Ytr = anomalies(Y[rows0], month[rows0], climatology(Y, month, rows0))
    Wset = Wm.build_all(Ytr, DATA, P["region_ids"]); names = list(Wset)
    W_unif = Wm.row_norm(np.ones((N, N)))
    Wset5 = dict(Wset); Wset5["uniform"] = W_unif
    Ya = anomalies(Y, month, climatology(Y, month, np.arange(n)))
    fits = {k: G.fit(Ya, W, X) for k, W in Wset5.items()}
    ref_dir = G.impacts(fits[REF]["beta"], Wset[REF])["direct"]; ref_ind = G.impacts(fits[REF]["beta"], Wset[REF])["indirect"]
    log = []

    # ---- A. uniform W as fifth candidate + H. dynamic-part margin
    rowsA = []
    for which in ["direct", "indirect"]:
        for a, b in itertools.combinations(list(Wset5), 2):
            th, psi = G.theta_and_influence(fits[a], fits[b], which)
            pr = infl(fits[REF], Wset[REF], which)
            r = rel_row(th, psi, ref_dir if which == "direct" else ref_ind, pr, rng)
            if which == "direct":
                r.update({k: v for k, v in rel_row(th, psi, ref_dir - 1.0, pr, rng, "dyn_").items()
                          if k.startswith("dyn_")})
            r.update(functional=which, W_A=a, W_B=b)
            rowsA.append(r)
    A_df = pd.DataFrame(rowsA); A_df.to_csv(RES / "extra_uniformW_equivalence.csv", index=False)
    est5 = []
    for k, f in fits.items():
        imp = G.impacts(f["beta"], Wset5[k])
        est5.append(dict(W=k, phi11_mean=G.phi11_mean(f["beta"]), direct=imp["direct"], indirect=imp["indirect"], rho_A=imp["rho"]))
    pd.DataFrame(est5).to_csv(RES / "extra_estimates_5W.csv", index=False)
    log.append(f"A/H done [{time.time()-t0:.0f}s]")

    # ---- B. common-factor augmentation (lagged cross-sectional mean)
    rowsB, estB = [], []
    fitsB = {k: fit_general(Ya, [W, W_unif], X) for k, W in Wset.items()}
    cols = [2, 3]   # columns of (W Y) and (mean Y) coefficients
    for k, f in fitsB.items():
        imp = impacts_multi(f["beta"], [Wset[k], W_unif], cols)
        estB.append(dict(W=k, phi11_mean=f["beta"][:, 2].mean(), phi_mean_coef=f["beta"][:, 3].mean(),
                         direct=imp["direct"], indirect=imp["indirect"], rho_A=imp["rho"],
                         indirect_W_only=impacts_multi(f["beta"], [Wset[k]], [2])["indirect"]))
    pd.DataFrame(estB).to_csv(RES / "extra_cce_estimates.csv", index=False)
    refB = pd.DataFrame(estB).set_index("W").loc[REF]
    for which in ["direct", "indirect"]:
        prB = G.influence_single(fitsB[REF], impact_grad_multi(fitsB[REF]["beta"], [Wset[REF], W_unif], cols, which))
        for a, b in itertools.combinations(names, 2):
            th, psi = theta_psi(fitsB[a], fitsB[b], [Wset[a], W_unif], [Wset[b], W_unif], cols, which)
            r = rel_row(th, psi, float(refB[which]), prB, rng)
            if which == "direct":
                r.update({k: v for k, v in rel_row(th, psi, float(refB[which]) - 1.0, prB, rng, "dyn_").items()
                          if k.startswith("dyn_")})
            r.update(functional=which, W_A=a, W_B=b); rowsB.append(r)
    pd.DataFrame(rowsB).to_csv(RES / "extra_cce_equivalence.csv", index=False)
    log.append(f"B done [{time.time()-t0:.0f}s]")

    # ---- C. similarity of spatial-lag series
    lags = {k: (Ya @ W.T) for k, W in Wset5.items()}
    sim = pd.DataFrame(index=list(Wset5), columns=list(Wset5), dtype=float)
    for a in Wset5:
        for b in Wset5:
            sim.loc[a, b] = float(np.mean([np.corrcoef(lags[a][:, i], lags[b][:, i])[0, 1] for i in range(N)]))
    sim.to_csv(RES / "extra_WY_similarity.csv")
    log.append(f"C done [{time.time()-t0:.0f}s]")

    # ---- D. J-test (Davidson-MacKinnon), both directions, per location HAC t and pooled
    rowsD = []
    for a, b in itertools.permutations(names, 2):
        fA, fB = fits[a], fits[b]
        Z = fA["Z"]; y = Ya[1:]; yhatB = fB["fitted"]           # (n-1, N)
        rej, tstats, scores, rho_i = 0, [], [], []
        num_pool = den_pool = 0.0
        for i in range(N):
            Zi = Z[:, i, :]; r = yhatB[:, i]
            # FWL: residualise y and r on Zi
            Qi = Zi.T @ Zi
            Py = y[:, i] - Zi @ np.linalg.solve(Qi, Zi.T @ y[:, i])
            Pr = r - Zi @ np.linalg.solve(Qi, Zi.T @ r)
            rho = float(Pr @ Py / (Pr @ Pr)); u = Py - rho * Pr
            sc = Pr * u                                           # score for rho
            se = np.sqrt(I.hac_lrv(sc) / sc.size) / (Pr @ Pr / sc.size)
            t = rho / se; tstats.append(t); rho_i.append(rho); rej += int(abs(t) > 1.96)
            num_pool += Pr @ Py; den_pool += Pr @ Pr; scores.append((Pr, Py))
        rho_pool = num_pool / den_pool
        sc_pool = np.sum([Pr * (Py - rho_pool * Pr) for Pr, Py in scores], axis=0)   # per-month summed score
        den = den_pool / sc_pool.size
        L = I.politis_white(sc_pool)
        se_pool = I.mbb_se(sc_pool, L, B, rng) / den
        rowsD.append(dict(maintained=a, rival=b, share_reject_5pct=rej / N, mean_rho=float(np.mean(rho_i)),
                          rho_pooled=rho_pool, se_pooled_mbb=se_pool, t_pooled=rho_pool / se_pool,
                          reject_pooled=bool(abs(rho_pool / se_pool) > 1.96)))
    pd.DataFrame(rowsD).to_csv(RES / "extra_jtest.csv", index=False)
    log.append(f"D done [{time.time()-t0:.0f}s]")

    # ---- E. benchmark RMSEs (rolling, m=240)
    targets = np.arange(M_WINDOW, n); sq_clim = np.empty((targets.size, idn.size)); sq_ar = np.empty_like(sq_clim)
    for j, t in enumerate(targets):
        rows = np.arange(t - M_WINDOW, t); clim = climatology(Y, month, rows)
        Yw = anomalies(Y[rows], month[rows], clim); Xw = X[rows]
        sq_clim[j] = (Y[t, idn] - clim[month[t] - 1][idn]) ** 2
        f = fit_general(Yw, [], Xw)                                # AR(1) + climate, no W
        y_prev = Yw[-1]; z = np.stack([np.ones(N), y_prev] + [np.repeat(X[t - 1, q], N) for q in range(X.shape[1])], axis=1)
        pred = np.einsum("ik,ik->i", z, f["beta"]) + clim[month[t] - 1]
        sq_ar[j] = (Y[t, idn] - pred[idn]) ** 2
    bench = [dict(model="climatology only", rmse_mm=float(np.sqrt(sq_clim.mean()))),
             dict(model="AR(1) + climate, no W", rmse_mm=float(np.sqrt(sq_ar.mean())))]
    fcx = {k: G.rolling_forecast(Y, month, W, X, m=M_WINDOW, score_idx=idn) for k, W in Wset5.items()}
    for k, v in fcx.items():
        bench.append(dict(model=f"GSTARX {k}", rmse_mm=v["rmse"]))
    # RMSE equivalence: GSTARX(W) vs AR(1) without W (does the spatial term matter at all?)
    for k, v in fcx.items():
        th, dt = I.rmse_diff_sequence(v["loss"], sq_ar.mean(axis=1)); r = eq_row(th, dt, 1.0, rng)
        bench.append(dict(model=f"GSTARX {k} minus AR(1): theta", rmse_mm=th, ci90_lo=r["ci90_lo"], ci90_hi=r["ci90_hi"]))
    pd.DataFrame(bench).to_csv(RES / "extra_benchmarks.csv", index=False)
    log.append(f"E done [{time.time()-t0:.0f}s]")

    # ---- F. sensitivity: KNN k and window m
    cent = pd.read_csv(DATA / "regions.csv", float_precision="round_trip").set_index("region_id")
    lat = np.array([cent.loc[r, "lat"] for r in P["region_ids"]]); lon = np.array([cent.loc[r, "lon"] for r in P["region_ids"]])
    dist = Wm.haversine_km(lat, lon); np.fill_diagonal(dist, np.inf)
    rowsF = []
    for k in (3, 5, 8):
        Wk = np.zeros((N, N))
        for i in range(N):
            nb = np.argsort(dist[i])[:k]
            Wk[i, nb] = 1.0 / dist[i, nb]        # inverse-distance weights within the k neighbours,
        Wk = Wm.row_norm(Wk)                      # as in knn_weights.csv (k=5 reproduces the main matrix)
        fk = G.fit(Ya, Wk, X); fk["W"] = Wk; imp = G.impacts(fk["beta"], Wk)
        th, psi = G.theta_and_influence(fk, fits["training-correlation"], "indirect")
        r = rel_row(th, psi, imp["indirect"], infl(fk, Wk, "indirect"), rng)
        rowsF.append(dict(setting=f"KNN k={k} (centroid-based)", direct=imp["direct"], indirect=imp["indirect"], rho_A=imp["rho"],
                          theta_vs_corr_indirect=th, c_min=r["c_min"],
                          overlap_with_a_geo=float(((Wk > 0) & (Wset[REF] > 0)).sum() / (Wset[REF] > 0).sum())))
    for m in (180, 240, 300):
        fa = G.rolling_forecast(Y, month, Wset[REF], X, m=m, score_idx=idn)
        fb = G.rolling_forecast(Y, month, Wset["training-correlation"], X, m=m, score_idx=idn)
        th, dt = I.rmse_diff_sequence(fa["loss"], fb["loss"]); r = eq_row(th, dt, 1.0, rng)
        rowsF.append(dict(setting=f"window m={m}: KNN vs corr RMSE", rmse_knn=fa["rmse"], rmse_corr=fb["rmse"], M=fa["loss"].size,
                          theta_rmse=th, min_delta=r["min_delta"], equivalent_1mm=r["equivalent"]))
    pd.DataFrame(rowsF).to_csv(RES / "extra_sensitivity.csv", index=False)
    log.append(f"F done [{time.time()-t0:.0f}s]")

    # ---- G. in-sample fits on months 241-484 only
    sub = np.arange(M_WINDOW, n); Ysub = anomalies(Y[sub], month[sub], climatology(Y, month, sub)); Xsub = X[sub]
    fitsG = {k: G.fit(Ysub, W, Xsub) for k, W in Wset.items()}
    rowsG = []
    impG = {k: G.impacts(f["beta"], Wset[k]) for k, f in fitsG.items()}
    for k in fitsG:
        rowsG.append(dict(W=k, kind="estimate", direct=impG[k]["direct"], indirect=impG[k]["indirect"], rho_A=impG[k]["rho"]))
    for which in ["direct", "indirect"]:
        pr = infl(fitsG[REF], Wset[REF], which); tau = impG[REF][which]
        for a, b in itertools.combinations(names, 2):
            th, psi = G.theta_and_influence(fitsG[a], fitsG[b], which)
            r = rel_row(th, psi, tau, pr, rng)
            row = dict(W=f"{a} vs {b}", kind=which, theta_hat=th, se_mbb=r["se_mbb"], c_min=r["c_min"],
                       equivalent=r["equivalent"], equivalent_hac=r["equivalent_hac"])
            if which == "direct":
                rd = rel_row(th, psi, tau - 1.0, pr, rng, "dyn_")
                row.update(dyn_c_min=rd["dyn_c_min"], dyn_equivalent=rd["dyn_equivalent"])
            rowsG.append(row)
    pd.DataFrame(rowsG).to_csv(RES / "extra_subsample_241_484.csv", index=False)
    log.append(f"G done [{time.time()-t0:.0f}s]")

    # ---- I. standardised anomalies: divide by the calendar-month standard deviation
    #      (full sample), which removes the seasonal pattern in the variance
    sd_c = np.array([[Ya[month == c, i].std(ddof=1) for i in range(N)] for c in range(1, 13)])
    Ys = Ya / sd_c[month - 1]
    fitsI = {k: G.fit(Ys, W, X) for k, W in Wset.items()}
    impI = {k: G.impacts(f["beta"], Wset[k]) for k, f in fitsI.items()}
    rowsI = [dict(W=k, kind="estimate", direct=impI[k]["direct"], indirect=impI[k]["indirect"], rho_A=impI[k]["rho"],
                  phi11_mean=G.phi11_mean(fitsI[k]["beta"])) for k in fitsI]
    for which in ["direct", "indirect"]:
        pr = infl(fitsI[REF], Wset[REF], which); tau = impI[REF][which]
        for a, b in itertools.combinations(names, 2):
            th, psi = G.theta_and_influence(fitsI[a], fitsI[b], which)
            r = rel_row(th, psi, tau, pr, rng)
            row = dict(W=f"{a} vs {b}", kind=which, theta_hat=th, se_mbb=r["se_mbb"], c_min=r["c_min"],
                       equivalent=r["equivalent"], equivalent_hac=r["equivalent_hac"])
            if which == "direct":
                rd = rel_row(th, psi, tau - 1.0, pr, rng, "dyn_")
                row.update(dyn_c_min=rd["dyn_c_min"], dyn_equivalent=rd["dyn_equivalent"])
            rowsI.append(row)
    pd.DataFrame(rowsI).to_csv(RES / "extra_standardised.csv", index=False)
    log.append(f"I done [{time.time()-t0:.0f}s]")

    with open(RES / "extra_summary.md", "w", encoding="utf-8") as fh:
        fh.write(f"# Extra analyses (generated {time.strftime('%Y-%m-%d %H:%M')})\n\n")
        fh.write("## A/H. Five candidates incl. uniform W; dynamic-part margin\n\n" + pd.DataFrame(est5).round(3).to_markdown(index=False) + "\n\n")
        fh.write(A_df.round(3).to_markdown(index=False) + "\n\n")
        fh.write("## B. Common-factor augmentation (lagged cross-sectional mean)\n\n" + pd.DataFrame(estB).round(3).to_markdown(index=False) + "\n\n" + pd.DataFrame(rowsB).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## C. Mean correlation between spatial-lag series\n\n" + sim.round(3).to_markdown() + "\n\n")
        fh.write("## D. J-test\n\n" + pd.DataFrame(rowsD).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## E. Benchmarks\n\n" + pd.DataFrame(bench).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## F. Sensitivity\n\n" + pd.DataFrame(rowsF).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## G. Months 241-484 only\n\n" + pd.DataFrame(rowsG).round(3).to_markdown(index=False) + "\n\n")
        fh.write("## I. Standardised anomalies\n\n" + pd.DataFrame(rowsI).round(3).to_markdown(index=False) + "\n\n")
        fh.write("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
