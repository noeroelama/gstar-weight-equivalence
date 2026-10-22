"""Simulation study for the equivalence test (Section 5 of the paper).

Experiment A (size and power of the TOST decision).
  DGP: heterogeneous GSTAR(1;1), N=18, true W0 = geographic KNN(k=3), spatially
  correlated Gaussian innovations (exponential, range 0.6). Candidate W_A = W0.
  Candidate W_B(lam) = row-normalised (1-lam) W0 + lam W_wrong, with lam chosen
  (root-finding on the closed-form population theta) so that |theta|/delta hits
  the grid {0, .25, .5, .75, .9, 1, 1.1, 1.25, 1.5, 2}. Quantities: phi11_mean and
  the long-run average direct effect. T in {120, 240, 480}. R replications; B
  bootstrap draws. For each replication we record whether equivalence is declared
  at alpha=0.05 with each standard-error estimator, and whether the two-sided 95%
  interval covers theta.
Experiment E (relative margin, as in the application): Part A design, margin
  c * tau(W_A) with tau(W_A) estimated; inference.tost_relative.
Experiment B (coverage across spatial-correlation strength x N).
  Same DGP family, lam fixed at 0.6, corr_range in {0.3, 0.6, 0.9}, N in {12, 24},
  T in {120, 240, 480}; coverage and mean width of the 95% interval for theta.

Standard-error estimators compared:
  mbb_pw   : vector MBB on the influence sequence, Politis-White block length
  mbb_L4   : vector MBB, fixed L=4
  mbb_L12  : vector MBB, fixed L=12
  hac      : Bartlett HAC with Newey-West bandwidth
  nodeind  : node-independent block bootstrap (breaks the cross-section), L=PW

Outputs: ../results/simA_size_power.csv, ../results/simB_coverage.csv, log.
"""
from __future__ import annotations
import sys, time, argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gstarx as G          # noqa: E402
import inference as I       # noqa: E402
import dgp as D             # noqa: E402

RES = HERE.parent / "results"
RES.mkdir(exist_ok=True)
ALPHA = 0.05
Z1 = norm.ppf(1 - ALPHA)          # one-sided for TOST
Z2 = norm.ppf(0.975)              # two-sided 95% CI
RATIOS = [0.0, 0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0]
TS = [120, 240, 480]
FUNCS = ["phi11_mean", "direct", "indirect"]


def pop_theta(dg, WA, WB, name):
    h = G.FUNCTIONALS[name][0]
    return h(dict(beta=D.pseudo_true(dg, WA), W=WA)) - h(dict(beta=D.pseudo_true(dg, WB), W=WB))


def lam_for_ratio(dg, WA, Wwrong, name, delta, ratio):
    f = lambda lam: abs(pop_theta(dg, WA, D.blend(WA, Wwrong, lam), name)) - ratio * delta
    if ratio == 0:
        return 0.0
    return brentq(f, 1e-6, 1.0, xtol=1e-6)


def se_all(fA, fB, name, psi, rng, B):
    """All standard errors for theta_hat = mean(psi)-based linearisation."""
    n = psi.size
    Lpw = I.politis_white(psi)
    out = dict(L_pw=Lpw)
    out["mbb_pw"] = I.mbb_se(psi, Lpw, B, rng)
    out["mbb_L4"] = I.mbb_se(psi, 4, B, rng)
    out["mbb_L12"] = I.mbb_se(psi, 12, B, rng)
    out["hac"] = I.hac_se(psi)
    gf = G.FUNCTIONALS[name][1]
    cA = np.einsum("tik,ik->ti", fA["g"], np.linalg.solve(fA["Q"], gf(fA)[..., None])[..., 0])
    cB = np.einsum("tik,ik->ti", fB["g"], np.linalg.solve(fB["Q"], gf(fB)[..., None])[..., 0])
    out["nodeind"] = I.node_indep_se(cA - cB, Lpw, B, rng)
    return out


METHODS = ["mbb_pw", "mbb_L4", "mbb_L12", "hac", "nodeind"]


Z025 = norm.ppf(1 - 0.025)


def one_rep(dg, WA, WB, name, T, theta0, delta, rng, B):
    sim = D.simulate(dg, T, rng)
    Y, X = (sim if isinstance(sim, tuple) else (sim, None))
    fA, fB = G.fit(Y, WA, X), G.fit(Y, WB, X)
    th, psi = G.theta_and_influence(fA, fB, name)
    ses = se_all(fA, fB, name, psi, rng, B)
    rec = dict(theta_hat=th, L_pw=ses["L_pw"], L_used=I.block_length(ses["L_pw"]),
               rho_max=max(G.spectral_radius(G.transition(fA["beta"], WA)),
                           G.spectral_radius(G.transition(fB["beta"], WB))))
    for mth in METHODS:
        se = ses[mth]
        lo, hi = th - Z1 * se, th + Z1 * se
        rec[f"eq_{mth}"] = int(lo > -delta and hi < delta)
        rec[f"cov_{mth}"] = int(abs(th - theta0) <= Z2 * se)
        rec[f"se_{mth}"] = se
        if mth in ("mbb_pw", "hac"):
            rec[f"eq025_{mth}"] = int(th - Z025 * se > -delta and th + Z025 * se < delta)
    return rec


def summarise(df):
    """Columns added to every results row from the per-replication records."""
    th = df.theta_hat.to_numpy()
    mad = np.median(np.abs(th - np.median(th)))
    out = dict(mean_theta_hat=th.mean(), sd_theta_hat=th.std(ddof=1), sd_robust=1.4826 * mad,
               mean_L_pw=df.L_pw.mean(), share_L2=float((df.L_used == 2).mean()),
               share_rho_ge_095=float((df.rho_max >= 0.95).mean()))
    for mth in METHODS:
        out[f"eq_{mth}"] = df[f"eq_{mth}"].mean()
        out[f"cov_{mth}"] = df[f"cov_{mth}"].mean()
        out[f"mean_se_{mth}"] = df[f"se_{mth}"].mean()
    for mth in ("mbb_pw", "hac"):
        out[f"eq025_{mth}"] = df[f"eq025_{mth}"].mean()
    return out


def experiment_A(R, B, seed=0, log=print):
    dg = D.make_dgp(N=18, k=3, corr_range=0.6, seed=seed)
    WA = dg["W0"]
    Wwrong = D.knn_from_coords(np.random.default_rng(seed + 101).random((18, 2)), 3)
    outA = RES / "simA_size_power.csv"
    rows = pd.read_csv(outA).to_dict("records") if outA.exists() else []
    done = {(r["functional"], round(r["ratio"], 4), int(r["T"])) for r in rows}
    for name in FUNCS:
        th_far = abs(pop_theta(dg, WA, D.blend(WA, Wwrong, 1.0), name))
        delta = round(th_far / 2.2, 4)          # so that ratio 2 is reachable
        log(f"[A] {name}: |theta(lam=1)| = {th_far:.4f}; delta = {delta}")
        for ratio in RATIOS:
            lam = lam_for_ratio(dg, WA, Wwrong, name, delta, ratio)
            WB = D.blend(WA, Wwrong, lam)
            theta0 = pop_theta(dg, WA, WB, name)
            for T in TS:
                if (name, round(ratio, 4), T) in done:
                    continue
                t0 = time.time()
                recs = []
                for r in range(R):
                    rng = np.random.default_rng(seed * 100003 + r)
                    recs.append(one_rep(dg, WA, WB, name, T, theta0, delta, rng, B))
                df = pd.DataFrame(recs)
                row = dict(functional=name, ratio=ratio, lam=lam, theta0=theta0, delta=delta,
                           T=T, R=R, B=B)
                row.update(summarise(df))
                rows.append(row)
                log(f"[A] {name} ratio={ratio:4.2f} T={T:3d}: eq(mbb_pw)={row['eq_mbb_pw']:.3f} "
                    f"eq(hac)={row['eq_hac']:.3f} eq(nodeind)={row['eq_nodeind']:.3f} "
                    f"cov(mbb_pw)={row['cov_mbb_pw']:.3f} sd={row['sd_theta_hat']:.4f} "
                    f"se={row['mean_se_mbb_pw']:.4f} [{time.time()-t0:.0f}s]")
                pd.DataFrame(rows).to_csv(RES / "simA_size_power.csv", index=False)
    return pd.DataFrame(rows)


def experiment_B(R, B, seed=0, log=print):
    outB = RES / "simB_coverage.csv"
    rows = pd.read_csv(outB).to_dict("records") if outB.exists() else []
    done = {(r["functional"], float(r["corr_range"]), int(r["N"]), int(r["T"])) for r in rows}
    for name in FUNCS:
        for corr_range in [0.3, 0.6, 0.9]:
            for N in [12, 24]:
                dg = D.make_dgp(N=N, k=3, corr_range=corr_range, seed=seed + N)
                WA = dg["W0"]
                Wwrong = D.knn_from_coords(np.random.default_rng(seed + 202 + N).random((N, 2)), 3)
                WB = D.blend(WA, Wwrong, 0.6)
                theta0 = pop_theta(dg, WA, WB, name)
                for T in TS:
                    if (name, float(corr_range), N, T) in done:
                        continue
                    t0 = time.time()
                    recs = [one_rep(dg, WA, WB, name, T, theta0, 1e9, np.random.default_rng(seed * 7919 + r), B)
                            for r in range(R)]
                    df = pd.DataFrame(recs)
                    row = dict(functional=name, corr_range=corr_range, N=N, T=T, R=R, B=B,
                               theta0=theta0)
                    row.update(summarise(df))
                    for mth in METHODS:
                        row[f"width_{mth}"] = (2 * Z2 * df[f"se_{mth}"]).mean()
                    rows.append(row)
                    log(f"[B] {name} corr={corr_range} N={N} T={T}: cov mbb_pw={row['cov_mbb_pw']:.3f} "
                        f"hac={row['cov_hac']:.3f} nodeind={row['cov_nodeind']:.3f} "
                        f"L12={row['cov_mbb_L12']:.3f} [{time.time()-t0:.0f}s]")
                    pd.DataFrame(rows).to_csv(RES / "simB_coverage.csv", index=False)
    return pd.DataFrame(rows)


def _ratio_sweep(dg, WA, Wwrong, ratios, Ts, funcs, R, B, seed, out_name, tag, log, extra=None):
    """Shared loop for Parts C and D: lambda chosen so that |theta|/delta hits each ratio."""
    out = RES / out_name
    rows = pd.read_csv(out).to_dict("records") if out.exists() else []
    done = {(r["functional"], round(r["ratio"], 4), int(r["T"])) for r in rows}
    for name in funcs:
        th_far = abs(pop_theta(dg, WA, D.blend(WA, Wwrong, 1.0), name))
        delta = round(th_far / 2.2, 4)
        log(f"[{tag}] {name}: |theta(lam=1)| = {th_far:.4f}; delta = {delta}")
        for ratio in ratios:
            lam = lam_for_ratio(dg, WA, Wwrong, name, delta, ratio)
            WB = D.blend(WA, Wwrong, lam)
            theta0 = pop_theta(dg, WA, WB, name)
            for T in Ts:
                if (name, round(ratio, 4), T) in done:
                    continue
                t0 = time.time()
                recs = [one_rep(dg, WA, WB, name, T, theta0, delta, np.random.default_rng(seed * 100003 + r), B)
                        for r in range(R)]
                row = dict(functional=name, ratio=ratio, lam=lam, theta0=theta0, delta=delta, T=T, R=R, B=B)
                row.update(extra or {})
                row.update(summarise(pd.DataFrame(recs)))
                rows.append(row)
                log(f"[{tag}] {name} ratio={ratio:4.2f} T={T:3d}: eq(mbb_pw)={row['eq_mbb_pw']:.3f} "
                    f"eq(hac)={row['eq_hac']:.3f} cov(mbb_pw)={row['cov_mbb_pw']:.3f} "
                    f"sd={row['sd_theta_hat']:.4f} se={row['mean_se_mbb_pw']:.4f} [{time.time()-t0:.0f}s]")
                pd.DataFrame(rows).to_csv(out, index=False)
    return pd.DataFrame(rows)


def experiment_C(R, B, seed=0, log=print):
    """Part C: neither candidate is the true matrix. The data come from a distance-band
    graph W0; W_A is the KNN(3) graph on the same locations; W_B(lambda) blends W_A towards
    an unrelated KNN(3) graph. Same N, coefficients and innovations as Part A."""
    dg = D.make_dgp(N=18, k=3, corr_range=0.6, seed=seed, graph="band")
    WA = D.knn_from_coords(dg["coords"], 3)
    Wwrong = D.knn_from_coords(np.random.default_rng(seed + 101).random((18, 2)), 3)
    log(f"[C] rho(A0) = {dg['rho']:.3f}; W_A shares {int(((WA > 0) & (dg['W0'] > 0)).sum())} of "
        f"{int((WA > 0).sum())} links with W0")
    return _ratio_sweep(dg, WA, Wwrong, [0.25, 0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5], TS, FUNCS,
                        R, B, seed, "simC_both_misspecified.csv", "C", log)


def experiment_D(R, B, seed=0, log=print):
    """Part D: application-sized design. N = 40 locations, two exogenous AR(1) variables
    entering with one lag, Student-t(5) innovations scaled to the same covariance, data
    from a distance-band graph, W_A = KNN(5) (misspecified), W_B(lambda) blends W_A towards
    an unrelated KNN(5) graph."""
    dg = D.make_dgp(N=40, k=5, corr_range=0.3, seed=seed + 40, graph="band", q=2, innov="t", df=5)
    WA = D.knn_from_coords(dg["coords"], 5)
    Wwrong = D.knn_from_coords(np.random.default_rng(seed + 140).random((40, 2)), 5)
    log(f"[D] rho(A0) = {dg['rho']:.3f}")
    return _ratio_sweep(dg, WA, Wwrong, [0.5, 1.0, 1.5], [240, 480], FUNCS,
                        R, B, seed, "simD_application_size.csv", "D", log,
                        extra=dict(N=40, q=2, innov="t5"))


def experiment_E(R, B, seed=0, log=print):
    """Part E: the relative margin used in the application. Same DGP, candidates and
    lambda values as Part A; the margin is c * tau(W_A) with c = delta / tau_pop(W_A), so
    that the population boundary |theta| = c tau(W_A) is the boundary of Part A. Both
    theta_hat and tau_hat(W_A) are estimated, and the test is inference.tost_relative."""
    dg = D.make_dgp(N=18, k=3, corr_range=0.6, seed=seed)
    WA = dg["W0"]
    Wwrong = D.knn_from_coords(np.random.default_rng(seed + 101).random((18, 2)), 3)
    out = RES / "simE_relative_margin.csv"
    rows = pd.read_csv(out).to_dict("records") if out.exists() else []
    done = {(r["functional"], round(r["ratio"], 4), int(r["T"])) for r in rows}
    for name in FUNCS:
        h, gradf = G.FUNCTIONALS[name]
        th_far = abs(pop_theta(dg, WA, D.blend(WA, Wwrong, 1.0), name))
        delta = round(th_far / 2.2, 4)
        tau_pop = h(dict(beta=D.pseudo_true(dg, WA), W=WA))
        c = delta / abs(tau_pop)
        log(f"[E] {name}: tau_pop(W_A) = {tau_pop:.4f}; delta = {delta}; c = {c:.4f}")
        for ratio in [0.5, 1.0, 1.5]:
            lam = lam_for_ratio(dg, WA, Wwrong, name, delta, ratio)
            WB = D.blend(WA, Wwrong, lam)
            theta0 = pop_theta(dg, WA, WB, name)
            for T in [240, 480]:
                if (name, round(ratio, 4), T) in done:
                    continue
                t0 = time.time()
                recs = []
                for r in range(R):
                    rng = np.random.default_rng(seed * 100003 + r)
                    sim = D.simulate(dg, T, rng)
                    Y, X = (sim if isinstance(sim, tuple) else (sim, None))
                    fA, fB = G.fit(Y, WA, X), G.fit(Y, WB, X)
                    th, psi = G.theta_and_influence(fA, fB, name)
                    tau_hat = h(fA)
                    psi_ref = G.influence_single(fA, gradf(fA))
                    sgn = np.sign(tau_hat) if tau_hat != 0 else 1.0
                    m = I.tost_relative(th, psi, sgn * tau_hat, sgn * psi_ref, c=c, B=B, rng=rng)
                    hv = I.tost_relative(th, psi, sgn * tau_hat, sgn * psi_ref, c=c,
                                         V=I.hac_cov2(psi, sgn * psi_ref))
                    recs.append(dict(eq_rel_mbb=int(m["equivalent_050"]), eq025_rel_mbb=int(m["equivalent_025"]),
                                     eq_rel_hac=int(hv["equivalent_050"]), eq025_rel_hac=int(hv["equivalent_025"]),
                                     L_used=m["L_used"]))
                df = pd.DataFrame(recs)
                row = dict(functional=name, ratio=ratio, lam=lam, theta0=theta0, delta=delta,
                           tau_pop=tau_pop, c=c, T=T, R=R, B=B)
                row.update({k: df[k].mean() for k in ["eq_rel_mbb", "eq025_rel_mbb", "eq_rel_hac", "eq025_rel_hac"]})
                row["share_L2"] = float((df.L_used == 2).mean())
                rows.append(row)
                log(f"[E] {name} ratio={ratio:4.2f} T={T:3d}: eq_rel(mbb)={row['eq_rel_mbb']:.3f} "
                    f"eq_rel(hac)={row['eq_rel_hac']:.3f} a025(mbb)={row['eq025_rel_mbb']:.3f} [{time.time()-t0:.0f}s]")
                pd.DataFrame(rows).to_csv(out, index=False)
    return pd.DataFrame(rows)


def write_design(seed=0):
    """Design constants quoted in Section 5, written to a results file."""
    lines = []
    for tag, kw in [("A", dict(N=18, k=3, corr_range=0.6, seed=seed)),
                    ("C", dict(N=18, k=3, corr_range=0.6, seed=seed, graph="band")),
                    ("D", dict(N=40, k=5, corr_range=0.3, seed=seed + 40, graph="band", q=2, innov="t"))]:
        dg = D.make_dgp(**kw)
        lines.append(f"Part {tag}: N={dg['N']} rho(A0)={dg['rho']:.3f} phi10 in [{dg['phi10'].min():.3f},"
                     f"{dg['phi10'].max():.3f}] phi11 in [{dg['phi11'].min():.3f},{dg['phi11'].max():.3f}] "
                     f"links in W0={int((dg['W0'] > 0).sum())}")
    dg = D.make_dgp(N=18, k=3, corr_range=0.6, seed=seed)
    Wwrong = D.knn_from_coords(np.random.default_rng(seed + 101).random((18, 2)), 3)
    lines.append(f"Part A: W_wrong shares {int(((Wwrong > 0) & (dg['W0'] > 0)).sum())} of "
                 f"{int((dg['W0'] > 0).sum())} links with W0")
    for N in (12, 24):
        for r in (0.3, 0.6, 0.9):
            dg = D.make_dgp(N=N, k=3, corr_range=r, seed=seed + N)
            lines.append(f"Part B: N={N} r={r} rho(A0)={dg['rho']:.3f}")
    (RES / "sim_design.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", default="AB")
    ap.add_argument("--R", type=int, default=1000)
    ap.add_argument("--B", type=int, default=399)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    logf = open(RES / f"sim_log_{a.exp}.txt", "a")

    def log(s):
        print(s, flush=True)
        logf.write(s + "\n"); logf.flush()

    log(f"=== sim_study start {time.strftime('%Y-%m-%d %H:%M:%S')} exp={a.exp} R={a.R} B={a.B}")
    if "A" in a.exp:
        experiment_A(a.R, a.B, a.seed, log)
    if "B" in a.exp:
        experiment_B(a.R, a.B, a.seed, log)
    if "C" in a.exp:
        experiment_C(a.R, a.B, a.seed, log)
    if "D" in a.exp:
        experiment_D(a.R, a.B, a.seed, log)
    if "E" in a.exp:
        experiment_E(a.R, a.B, a.seed, log)
    write_design(a.seed)
    log(f"=== done {time.strftime('%Y-%m-%d %H:%M:%S')}")
