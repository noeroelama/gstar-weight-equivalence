"""Sanity checks for the core functions. Run: python test_core.py
 1. impact gradients vs finite differences
 2. pseudo-true coefficients (closed form) vs OLS on a very long simulated series
 3. se of theta_hat (MBB on influence sequence, HAC) vs Monte-Carlo sd
 4. application smoke test: anomalies, W set, in-sample fits, rolling forecast
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gstarx as G                       # noqa: E402
import inference as I                    # noqa: E402
import dgp as D                          # noqa: E402

rng = np.random.default_rng(1)
print("== 1. impact gradient vs finite differences")
dg = D.make_dgp(N=10, seed=3)
beta = np.column_stack([np.zeros(10), dg["phi10"], dg["phi11"]])
for which in ["direct", "indirect"]:
    Gan = G.impact_grad(beta, dg["W0"], which)
    Gfd = np.zeros_like(beta)
    h = 1e-6
    for i in range(10):
        for c in (1, 2):
            bp = beta.copy(); bp[i, c] += h
            bm = beta.copy(); bm[i, c] -= h
            Gfd[i, c] = (G.impacts(bp, dg["W0"])[which] - G.impacts(bm, dg["W0"])[which]) / (2 * h)
    print(f"  {which:9s} max |analytic - FD| = {np.abs(Gan - Gfd).max():.2e}")

print("== 2. pseudo-true (closed form) vs long-sample OLS under a WRONG W")
dg = D.make_dgp(N=12, seed=5)
Wwrong = D.knn_from_coords(np.random.default_rng(99).random((12, 2)), 3)
WB = D.blend(dg["W0"], Wwrong, 0.6)
Y = D.simulate(dg, 400_000, rng)
f = G.fit(Y, WB)
bstar = D.pseudo_true(dg, WB)
print(f"  max |beta_hat - beta*| over 12x3 = {np.abs(f['beta'] - bstar).max():.4f} "
      f"(sampling noise ~ {1/np.sqrt(400_000):.4f})")
print(f"  rho(A*(W_B)) = {G.spectral_radius(G.transition(bstar, WB)):.3f}, rho(A0) = {dg['rho']:.3f}")

print("== 3. se of theta_hat: MBB(influence) and HAC vs Monte-Carlo sd")
dg = D.make_dgp(N=12, seed=5, corr_range=0.6)
WA, WBm = dg["W0"], WB
T, R, B = 300, 400, 300
for name in ["phi11_mean", "direct"]:
    th = np.empty(R); se_mbb = np.empty(R); se_hac = np.empty(R); se_ni = np.empty(R)
    t0 = time.time()
    for r in range(R):
        rr = np.random.default_rng(1000 + r)
        Y = D.simulate(dg, T, rr)
        fA, fB = G.fit(Y, WA), G.fit(Y, WBm)
        th[r], psi = G.theta_and_influence(fA, fB, name)
        L = I.politis_white(psi)
        se_mbb[r] = I.mbb_se(psi, L, B, rr) / 1.0 * 1.0
        se_hac[r] = I.hac_se(psi)
        # per-location contributions for the node-independent comparator
        hA, gA = G.FUNCTIONALS[name]; hB, gB = G.FUNCTIONALS[name]
        cA = np.einsum("tik,ik->ti", fA["g"], np.linalg.solve(fA["Q"], gA(fA)[..., None])[..., 0])
        cB = np.einsum("tik,ik->ti", fB["g"], np.linalg.solve(fB["Q"], gB(fB)[..., None])[..., 0])
        se_ni[r] = I.node_indep_se(cA - cB, L, B, rr)
    n = fA["n"]
    print(f"  {name:10s} MC sd(theta_hat) = {th.std(ddof=1):.4f} | mean se: "
          f"MBB/sqrt(n) = {np.mean(se_mbb)/1:.4f}  HAC = {np.mean(se_hac):.4f}  "
          f"node-indep = {np.mean(se_ni):.4f}   [{time.time()-t0:.1f}s]")

print("== 4. application smoke test")
sys.path.insert(0, str(HERE))
from data_io import load_panel, climatology, anomalies   # noqa: E402
import weights as Wm                                       # noqa: E402
P = load_panel()
Y, X, month, idn = P["Y"], P["X"], P["month"], P["idn"]
print(f"  panel {Y.shape}, X {X.shape}, scored locations {idn.size}")
m = 240
Wset = Wm.build_all(anomalies(Y[:m], month[:m], climatology(Y, month, np.arange(m))),
                    HERE.parent / "data", P["region_ids"])
clim_full = climatology(Y, month, np.arange(Y.shape[0]))
Ya = anomalies(Y, month, clim_full)
for nm, W in Wset.items():
    f = G.fit(Ya, W, X)
    imp = G.impacts(f["beta"], W)
    print(f"  {nm:22s} phi11_mean={G.phi11_mean(f['beta']):+.3f}  dir={imp['direct']:.3f} "
          f"ind={imp['indirect']:.3f} rho(A)={imp['rho']:.3f}  #phi11<0: {(f['beta'][:,2]<0).sum()}")
t0 = time.time()
rA = G.rolling_forecast(Y, month, Wset["geographic KNN"], X, m=m, score_idx=idn)
rB = G.rolling_forecast(Y, month, Wset["training-correlation"], X, m=m, score_idx=idn)
th, dt = I.rmse_diff_sequence(rA["loss"], rB["loss"])
L = I.politis_white(dt)
res = I.tost(th, I.mbb_se(dt, L, 999, rng), delta=1.0)
print(f"  rolling m={m}: M={rA['loss'].size} months, RMSE KNN={rA['rmse']:.2f} corr={rB['rmse']:.2f} mm; "
      f"theta={th:+.3f} se={res['se']:.3f} L_PW={L:.1f} min_delta={res['min_delta']:.2f}  [{time.time()-t0:.1f}s]")
