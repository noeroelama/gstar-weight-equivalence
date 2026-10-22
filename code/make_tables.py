"""Build the LaTeX table bodies for Sections 5-6 and Appendix C from the CSVs in ../results.
Every number in the manuscript tables comes from here; nothing is typed by hand.
Output: ../manuscript/tables/*.tex (each file is a complete tabular environment) and
../manuscript/tables/sim_constants.tex (design constants quoted in Section 5).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
TAB = HERE.parent / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SHORT = {"geographic KNN": "KNN", "inverse-distance": "Inv.\\ dist.",
         "training-correlation": "Corr.", "graphical-lasso": "G-lasso", "uniform": "Uniform"}
QNAME = {"phi11_mean": r"$\bar\phi_{11}$", "direct": "direct effect", "indirect": "indirect effect"}
FUNCS = ["phi11_mean", "direct", "indirect"]


def f3(x):
    return f"{x:.3f}"


def f2(x):
    return f"{x:.2f}".replace("-0.00", "0.00")


def pct(x):
    return "--" if pd.isna(x) else f"{100 * x:.1f}"


def yn(x):
    return "yes" if bool(x) else "no"


def write(name, head, rows):
    body = head + "\n".join(rows) + "\n" + r"\bottomrule" + "\n" + r"\end{tabular}" + "\n"
    (TAB / name).write_text(body, encoding="utf-8")
    print("wrote", name)

# ---------------------------------------------------------------- simulation


def _cell(D, func, ratio, T):
    g = D[(D.functional == func) & np.isclose(D.ratio, ratio) & (D["T"] == T)]
    return g.iloc[0] if len(g) else None


def _boundary_reps():
    """Per-replication records at the boundary (sim_boundary_reps.py), if present."""
    files = sorted(RES.glob("sim_boundary_reps_*.csv"))
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True) if files else None


BREPS = None


def _spread(part, func, T, b):
    """MC sd (usual, robust) of theta-hat and median se (MBB, HAC) at the boundary."""
    global BREPS
    if BREPS is None:
        BREPS = _boundary_reps()
    if BREPS is not None:
        g = BREPS[(BREPS.part == part) & (BREPS.functional == func) & (BREPS["T"] == T)]
        if len(g):
            th = g.theta_hat.to_numpy()
            rob = 1.4826 * np.median(np.abs(th - np.median(th)))
            return th.std(ddof=1), rob, g.se_mbb_pw.median(), g.se_hac.median()
    return b.sd_theta_hat, b.sd_robust, b.mean_se_mbb_pw, b.mean_se_hac


def table_sim_summary(csv, name, Ts, power_ratio=0.5, beyond=1.25, nodeind=True, part="A"):
    """One row per quantity and T: power inside the region, rate at the boundary (several
    standard errors, alpha = 0.05 and 0.025), rate beyond it, coverage at the boundary, and
    the Monte-Carlo sd of theta-hat (usual and robust) against the median estimated standard
    errors (from the per-replication boundary records when available)."""
    D = pd.read_csv(RES / csv)
    rows = []
    for func in FUNCS:
        if func not in set(D.functional):
            continue
        first = True
        for T in Ts:
            p, b, o = _cell(D, func, power_ratio, T), _cell(D, func, 1.0, T), _cell(D, func, beyond, T)
            if b is None:
                continue
            cells = [QNAME[func] if first else "", str(int(T)), f3(p.eq_mbb_pw), f3(b.eq_mbb_pw), f3(b.eq_hac)]
            if nodeind:
                cells.append(f3(b.eq_nodeind))
            sd, rob, se_m, se_h = _spread(part, func, T, b)
            cells += [f3(b.eq025_mbb_pw), f3(o.eq_mbb_pw) if o is not None else "--", f3(b.cov_mbb_pw),
                      f"{sd:.4f}", f"{rob:.4f}", f"{se_m:.4f}", f"{se_h:.4f}"]
            rows.append(" & ".join(cells) + r" \\")
            first = False
        rows.append(r"\addlinespace")
    ncol_eq = 4 if nodeind else 3
    head = (r"\begin{tabular}{l r c " + "c" * ncol_eq + r" c c cc cc}" "\n" r"\toprule" "\n"
            r" & & power & \multicolumn{" + str(ncol_eq) + r"}{c}{rate at the boundary $|\theta|=\delta$} & rate & coverage & \multicolumn{2}{c}{MC sd of $\hat\theta$} & \multicolumn{2}{c}{median se} \\" "\n"
            r"\cmidrule(lr){4-" + str(3 + ncol_eq) + r"}\cmidrule(lr){" + str(6 + ncol_eq) + "-" + str(7 + ncol_eq) + r"}\cmidrule(lr){" + str(8 + ncol_eq) + "-" + str(9 + ncol_eq) + r"}" "\n"
            r"quantity & $T$ & $|\theta|/\delta=" + f"{power_ratio:g}" + r"$ & MBB & HAC & " + ("node-ind.\\ & " if nodeind else "")
            + r"MBB, $\alpha{=}0.025$ & $|\theta|/\delta=" + f"{beyond:g}" + r"$ & at boundary & usual & robust & MBB & HAC \\" "\n"
            r"\midrule" "\n")
    write(name, head, rows[:-1])


def table_simCD():
    """Parts C and D in one table: rows by part, quantity and T. The rate beyond the boundary
    is at |theta|/delta = 1.25 in Part C and 1.5 in Part D."""
    rows = []
    for part, csv, Ts, beyond in [("C", "simC_both_misspecified.csv", [120, 240, 480], 1.25),
                                  ("D", "simD_application_size.csv", [240, 480], 1.5)]:
        D = pd.read_csv(RES / csv)
        firstpart = True
        for func in FUNCS:
            first = True
            for T in Ts:
                p, b, o = _cell(D, func, 0.5, T), _cell(D, func, 1.0, T), _cell(D, func, beyond, T)
                if b is None:
                    continue
                sd, rob, se_m, se_h = _spread(part, func, T, b)
                cells = [part if firstpart else "", QNAME[func] if first else "", str(int(T)),
                         f3(p.eq_mbb_pw), f3(b.eq_mbb_pw), f3(b.eq_hac), f3(b.eq025_mbb_pw),
                         f3(o.eq_mbb_pw) if o is not None else "--", f3(b.cov_mbb_pw),
                         f"{sd:.4f}", f"{rob:.4f}", f"{se_m:.4f}", f"{se_h:.4f}"]
                rows.append(" & ".join(cells) + r" \\")
                first = firstpart = False
            rows.append(r"\addlinespace")
        rows[-1] = r"\midrule"
    head = (r"\begin{tabular}{l l r c ccc c c cc cc}" "\n" r"\toprule" "\n"
            r" & & & power & \multicolumn{3}{c}{rate at the boundary $|\theta|=\delta$} & rate & coverage & \multicolumn{2}{c}{MC sd of $\hat\theta$} & \multicolumn{2}{c}{median se} \\" "\n"
            r"\cmidrule(lr){5-7}\cmidrule(lr){10-11}\cmidrule(lr){12-13}" "\n"
            r"part & quantity & $T$ & $|\theta|/\delta=0.5$ & MBB & HAC & MBB, $\alpha{=}0.025$ & beyond & at boundary & usual & robust & MBB & HAC \\" "\n"
            r"\midrule" "\n")
    write("tab_simCD_summary.tex", head, rows[:-1])


def table_simA_full(func, label):
    """Appendix C: the whole ratio grid of Part A for one quantity."""
    A = pd.read_csv(RES / "simA_size_power.csv")
    A = A[(A.functional == func) & (A.ratio > 0)]
    rows = []
    for ratio, g in A.groupby("ratio"):
        g = g.sort_values("T")
        cells = [f"{ratio:.2f}"]
        for col in ["eq_mbb_pw", "eq_nodeind", "cov_mbb_pw"]:
            cells += [f3(v) for v in g[col]]
        rows.append(" & ".join(cells) + r" \\")
    head = (r"\begin{tabular}{l ccc ccc ccc}" "\n" r"\toprule" "\n"
            r" & \multicolumn{6}{c}{share declaring equivalence} & \multicolumn{3}{c}{coverage of the 95\% interval} \\" "\n"
            r"\cmidrule(lr){2-7}\cmidrule(lr){8-10}" "\n"
            r" & \multicolumn{3}{c}{vector MBB} & \multicolumn{3}{c}{node-independent} & \multicolumn{3}{c}{vector MBB} \\" "\n"
            r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}" "\n"
            r"$|\theta|/\delta$ & $T{=}120$ & 240 & 480 & 120 & 240 & 480 & 120 & 240 & 480 \\" "\n" r"\midrule" "\n")
    write(f"tab_simA_full_{label}.tex", head, rows)
    return float(A.delta.iloc[0])


def table_simB():
    B = pd.read_csv(RES / "simB_coverage.csv")
    rows = []
    for func in FUNCS:
        first = True
        for (cr, N), g in B[B.functional == func].groupby(["corr_range", "N"]):
            g = g.sort_values("T")
            cells = [QNAME[func] if first else "", f"{cr:.1f}", str(int(N)), f"{g.theta0.iloc[0]:+.3f}"]
            for mth in ["mbb_pw", "hac", "nodeind"]:
                cells += [f3(v) for v in g[f"cov_{mth}"]]
            rows.append(" & ".join(cells) + r" \\")
            first = False
        rows.append(r"\addlinespace")
    head = (r"\begin{tabular}{l cc r ccc ccc ccc}" "\n" r"\toprule" "\n"
            r" & & & & \multicolumn{3}{c}{vector MBB} & \multicolumn{3}{c}{HAC} & \multicolumn{3}{c}{node-independent} \\" "\n"
            r"\cmidrule(lr){5-7}\cmidrule(lr){8-10}\cmidrule(lr){11-13}" "\n"
            r"quantity & $r$ & $N$ & $\theta$ & $T{=}120$ & 240 & 480 & 120 & 240 & 480 & 120 & 240 & 480 \\" "\n"
            r"\midrule" "\n")
    write("tab_simB.tex", head, rows[:-1])


def table_simE():
    """Part E: relative margin c * tau_hat(W_A), next to the absolute margin of Part A."""
    E = pd.read_csv(RES / "simE_relative_margin.csv")
    A = pd.read_csv(RES / "simA_size_power.csv")
    rows = []
    for func in FUNCS:
        first = True
        for T in [240, 480]:
            p, b, o = (_cell(E, func, x, T) for x in (0.5, 1.0, 1.5))
            a = _cell(A, func, 1.0, T)
            if b is None:
                continue
            cells = [QNAME[func] if first else "", (f"{100 * b.c:.1f}" if first else ""), str(int(T)),
                     f3(p.eq_rel_mbb), f3(b.eq_rel_mbb), f3(b.eq_rel_hac), f3(b.eq025_rel_mbb),
                     f3(o.eq_rel_mbb), (f3(a.eq_mbb_pw) if a is not None else "--")]
            rows.append(" & ".join(cells) + r" \\")
            first = False
        rows.append(r"\addlinespace")
    head = (r"\begin{tabular}{l r r c ccc c c}" "\n" r"\toprule" "\n"
            r" & & & power & \multicolumn{3}{c}{rate at the boundary} & rate & absolute margin \\" "\n"
            r"\cmidrule(lr){5-7}" "\n"
            r"quantity & $c$ (\%) & $T$ & $|\theta|/\delta=0.5$ & MBB & HAC & MBB, $\alpha{=}0.025$ & $|\theta|/\delta=1.5$ & (Part A, MBB) \\" "\n"
            r"\midrule" "\n")
    write("tab_simE_summary.tex", head, rows[:-1])

# ---------------------------------------------------------------- application


def table_app_similarity():
    S = pd.read_csv(RES / "extra_WY_similarity.csv", index_col=0)
    names = list(S.index)
    rows = [" & ".join([SHORT[a]] + [f2(S.loc[a, b]) if j > i else "" for j, b in enumerate(names)]) + r" \\"
            for i, a in enumerate(names[:-1])]
    head = (r"\begin{tabular}{l" + "c" * len(names) + "}" "\n" r"\toprule" "\n" + " & "
            + " & ".join(SHORT[b] for b in names) + r" \\" "\n" r"\midrule" "\n")
    write("tab_app_similarity.tex", head, rows)


def table_app_estimates():
    E = pd.read_csv(RES / "app_insample_estimates.csv")
    rows = []
    for _, r in E.iterrows():
        rows.append(" & ".join([SHORT[r.W], f3(r.phi10_mean), f3(r.phi11_mean), str(int(r.n_phi11_negative)),
                                f"{r.direct:.3f} ({r.se_direct:.3f})", f"{r.indirect:.3f} ({r.se_indirect:.3f})",
                                f3(r.rho_A), f2(r.resid_sd_mm), f2(r.ljung_box_reject_share),
                                f2(r.mean_abs_resid_crosscorr)]) + r" \\")
    head = (r"\begin{tabular}{l ccc cc c ccc}" "\n" r"\toprule" "\n"
            r"$W$ & mean $\hat\phi_{10}$ & mean $\hat\phi_{11}$ & \#$\hat\phi_{11}{<}0$ & direct (se) & indirect (se) & $\rho(\hat A)$ & resid.\ sd & LB reject & $|$resid.\ corr.$|$ \\" "\n"
            r"\midrule" "\n")
    write("tab_app_estimates.tex", head, rows)


REL_HEAD = (r" & & & & & \multicolumn{2}{c}{relative to $\tau(\mathrm{KNN})$} & \multicolumn{2}{c}{relative to $\tau_{\mathrm{dir}}(\mathrm{KNN})-1$} \\" "\n"
            r"\cmidrule(lr){6-7}\cmidrule(lr){8-9}" "\n"
            r"quantity & pair & $\hat\theta$ & se & 90\% interval & $c_{\min}$ (\%) & equiv.\ at 10\%? & $c_{\min}$ (\%) & equiv.\ at 10\%? \\" "\n"
            r"\midrule" "\n")


def table_app_eq12():
    Q = pd.read_csv(RES / "app_equivalence_Q12.csv")
    rows = []
    for func in FUNCS:
        g = Q[Q.functional == func]
        first = True
        for _, r in g.iterrows():
            cells = [QNAME[func] if first else "", f"{SHORT[r.W_A]} vs {SHORT[r.W_B]}", f"{r.theta_hat:+.3f}", f3(r.se_mbb),
                     f"[{r.ci90_lo:+.3f}, {r.ci90_hi:+.3f}]", pct(r.c_min_ref), yn(r.equivalent_ref)]
            cells += [pct(r.c_min_dyn), yn(r.equivalent_dyn)] if func == "direct" else ["--", "--"]
            rows.append(" & ".join(cells) + r" \\")
            first = False
        rows.append(r"\addlinespace")
    write("tab_app_eq12.tex", r"\begin{tabular}{l l r r c r c r c}" "\n" r"\toprule" "\n" + REL_HEAD, rows[:-1])


def table_app_bench():
    """RMSE by matrix, with and without the climate indices, and the GSTARX difference from the
    AR(1) forecast with the indices but no spatial lag."""
    R = pd.read_csv(RES / "app_rmse.csv")
    P = R.pivot(index="W", columns="model", values="rmse_mm")
    Bm = pd.read_csv(RES / "extra_benchmarks.csv")
    rows = []
    for _, r in Bm[~Bm.model.str.contains("GSTARX")].iterrows():
        rows.append(f"{r.model} & -- & {f2(r.rmse_mm)} & -- & --" + r" \\")
    for w in list(SHORT):
        d = Bm[Bm.model == f"GSTARX {w} minus AR(1): theta"]
        x = Bm[Bm.model == f"GSTARX {w}"]
        if not len(x):
            continue
        gstar = f2(P.loc[w, "GSTAR"]) if w in P.index else "--"
        rows.append(" & ".join([f"GSTAR(X), {SHORT[w]}", gstar, f2(x.iloc[0].rmse_mm),
                                f"{d.iloc[0].rmse_mm:+.2f}", f"[{d.iloc[0].ci90_lo:+.2f}, {d.iloc[0].ci90_hi:+.2f}]"]) + r" \\")
    head = (r"\begin{tabular}{l c c r c}" "\n" r"\toprule" "\n"
            r" & \multicolumn{2}{c}{RMSE (mm)} & \multicolumn{2}{c}{with indices, minus AR(1)} \\" "\n"
            r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}" "\n"
            r"forecast & without indices & with indices & $\hat\theta$ (mm) & 90\% interval \\" "\n" r"\midrule" "\n")
    write("tab_app_bench.tex", head, rows)


def table_app_eq_rmse():
    Q = pd.read_csv(RES / "app_equivalence_rmse.csv")
    rows = []
    for model, lab in [("GSTARX", "with indices"), ("GSTAR", "without indices")]:
        g = Q[Q.model == model]
        first = True
        for _, r in g.iterrows():
            rows.append(" & ".join([lab if first else "", f"{SHORT[r.W_A]} vs {SHORT[r.W_B]}", f2(r.theta_rmse_mm),
                                    f2(r.se_mbb), f2(r.se_hac), f2(r.se_mbb_L12), f"[{r.ci90_lo:+.2f}, {r.ci90_hi:+.2f}]",
                                    yn(r.equivalent_mbb), f2(r.min_delta_mm), str(int(r.L_used))]) + r" \\")
            first = False
        rows.append(r"\addlinespace")
    head = (r"\begin{tabular}{l l r rrr c c r r}" "\n" r"\toprule" "\n"
            r"model & pair & $\hat\theta$ (mm) & se (MBB) & se (HAC) & se (MBB, $L{=}12$) & 90\% interval & equiv.\ at 1 mm? & $\delta_{\min}$ (mm) & $L$ \\" "\n"
            r"\midrule" "\n")
    write("tab_app_eq_rmse.tex", head, rows[:-1])


def table_app_cce():
    E = pd.read_csv(RES / "extra_cce_estimates.csv")
    F = pd.read_csv(RES / "extra_estimates_5W.csv").set_index("W")
    rows = [" & ".join([SHORT[r.W], f3(F.loc[r.W, "indirect"]), f3(r.phi11_mean), f3(r.phi_mean_coef), f3(r.indirect_W_only),
                        f3(r.indirect), f3(r.direct), f3(r.rho_A)]) + r" \\" for _, r in E.iterrows()]
    head = (r"\begin{tabular}{l c cccccc}" "\n" r"\toprule" "\n"
            r" & without mean term & \multicolumn{6}{c}{with the lagged mean of the other regions} \\" "\n"
            r"\cmidrule(lr){2-2}\cmidrule(lr){3-8}" "\n"
            r"$W$ & indirect & mean $\hat\phi_{11}$ & mean $\hat\phi_{m}$ & indirect ($W$ part) & indirect (total) & direct & $\rho(\hat A)$ \\" "\n"
            r"\midrule" "\n")
    write("tab_app_cce.tex", head, rows)
    Qc = pd.read_csv(RES / "extra_cce_equivalence.csv")
    rows = []
    for func in ["direct", "indirect"]:
        g = Qc[Qc.functional == func]
        first = True
        for _, r in g.iterrows():
            cells = [QNAME[func] if first else "", f"{SHORT[r.W_A]} vs {SHORT[r.W_B]}", f"{r.theta_hat:+.3f}", f3(r.se_mbb),
                     f"[{r.ci90_lo:+.3f}, {r.ci90_hi:+.3f}]", pct(r.c_min), yn(r.equivalent)]
            cells += [pct(r.dyn_c_min), yn(r.dyn_equivalent)] if func == "direct" else ["--", "--"]
            rows.append(" & ".join(cells) + r" \\")
            first = False
        rows.append(r"\addlinespace")
    write("tab_app_cce_eq.tex", r"\begin{tabular}{l l r r c r c r c}" "\n" r"\toprule" "\n" + REL_HEAD, rows[:-1])


def table_app_jtest():
    J = pd.read_csv(RES / "extra_jtest.csv")
    rows = [" & ".join([SHORT[r.maintained], SHORT[r.rival], f"{int(round(40 * r.share_reject_5pct))}/40", f"{r.rho_pooled:+.2f}", f2(r.se_pooled_mbb),
                        f"{r.t_pooled:+.1f}", yn(r.reject_pooled)]) + r" \\" for _, r in J.iterrows()]
    head = (r"\begin{tabular}{l l c r r r c}" "\n" r"\toprule" "\n"
            r"maintained $W$ & rival $W$ & locations rejecting & pooled $\hat\lambda$ & se & $t$ & rejects? \\" "\n" r"\midrule" "\n")
    write("tab_app_jtest.tex", head, rows)

# ---------------------------------------------------------------- Appendix C


def table_app_raw():
    E = pd.read_csv(RES / "app_insample_rawdata.csv")
    A = pd.read_csv(RES / "app_insample_estimates.csv").set_index("W")
    rows = []
    for _, r in E.iterrows():
        a = A.loc[r.W]
        rows.append(" & ".join([SHORT[r.W], f3(r.phi10_mean), f3(a.phi10_mean), f3(r.rho_A), f3(a.rho_A),
                                f3(r.indirect), f3(a.indirect), f2(r.ljung_box_reject_share), f2(a.ljung_box_reject_share)]) + r" \\")
    head = (r"\begin{tabular}{l cc cc cc cc}" "\n" r"\toprule" "\n"
            r" & \multicolumn{2}{c}{mean $\hat\phi_{10}$} & \multicolumn{2}{c}{$\rho(\hat A)$} & \multicolumn{2}{c}{indirect effect} & \multicolumn{2}{c}{LB reject share} \\" "\n"
            r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}\cmidrule(lr){6-7}\cmidrule(lr){8-9}" "\n"
            r"$W$ & raw & anomaly & raw & anomaly & raw & anomaly & raw & anomaly \\" "\n"
            r"\midrule" "\n")
    write("tab_app_raw.tex", head, rows)


def _eq_block(G, est_cols):
    """Estimates and pairwise tests for one robustness fit (subsample or standardised)."""
    est = G[G.kind == "estimate"]
    rows_e = [" & ".join([SHORT[r.W]] + [f3(r[c]) for c in est_cols]) + r" \\" for _, r in est.iterrows()]
    rows_q = []
    for func in ["direct", "indirect"]:
        g = G[G.kind == func]
        first = True
        for _, r in g.iterrows():
            a, b = r.W.split(" vs ")
            cells = [QNAME[func] if first else "", f"{SHORT[a]} vs {SHORT[b]}", f"{r.theta_hat:+.3f}", f3(r.se_mbb),
                     pct(r.c_min), yn(r.equivalent)]
            cells += [pct(r.dyn_c_min), yn(r.dyn_equivalent)] if func == "direct" else ["--", "--"]
            rows_q.append(" & ".join(cells) + r" \\")
            first = False
        rows_q.append(r"\addlinespace")
    return rows_e, rows_q[:-1]


def table_app_robust():
    Sv = pd.read_csv(RES / "extra_sensitivity.csv")
    rows = [" & ".join([r.setting.replace(" (centroid-based)", "").replace("KNN ", ""), f3(r.direct), f3(r.indirect), f3(r.rho_A),
                        f"{r.theta_vs_corr_indirect:+.3f}", pct(r.c_min)]) + r" \\" for _, r in Sv[Sv.setting.str.startswith("KNN")].iterrows()]
    head = (r"\begin{tabular}{l ccc r r}" "\n" r"\toprule" "\n"
            r"neighbours & direct & indirect & $\rho(\hat A)$ & $\hat\theta$, indirect vs Corr.\ & $c_{\min}$ (\%) \\" "\n" r"\midrule" "\n")
    write("tab_app_sens_k.tex", head, rows)
    rows = []
    for _, r in Sv[Sv.setting.str.startswith("window")].iterrows():
        m = r.setting.split("=")[1].split(":")[0]
        rows.append(" & ".join([m, str(int(r.M)), f2(r.rmse_knn), f2(r.rmse_corr), f"{r.theta_rmse:+.2f}", f2(r.min_delta), yn(r.equivalent_1mm)]) + r" \\")
    head = (r"\begin{tabular}{c c cc r r c}" "\n" r"\toprule" "\n"
            r"$m$ & $M$ & RMSE KNN & RMSE Corr.\ & $\hat\theta$ (mm) & $\delta_{\min}$ (mm) & equiv.\ at 1 mm? \\" "\n" r"\midrule" "\n")
    write("tab_app_sens_m.tex", head, rows)
    head_e = (r"\begin{tabular}{l ccc}" "\n" r"\toprule" "\n" r"$W$ & direct & indirect & $\rho(\hat A)$ \\" "\n" r"\midrule" "\n")
    head_q = (r"\begin{tabular}{l l r r r c r c}" "\n" r"\toprule" "\n"
              r" & & & & \multicolumn{2}{c}{relative to $\tau(\mathrm{KNN})$} & \multicolumn{2}{c}{relative to $\tau_{\mathrm{dir}}(\mathrm{KNN})-1$} \\" "\n"
              r"\cmidrule(lr){5-6}\cmidrule(lr){7-8}" "\n"
              r"quantity & pair & $\hat\theta$ & se & $c_{\min}$ (\%) & equiv.? & $c_{\min}$ (\%) & equiv.? \\" "\n" r"\midrule" "\n")
    for csv, tag in [("extra_subsample_241_484.csv", "sub"), ("extra_standardised.csv", "std")]:
        e, q = _eq_block(pd.read_csv(RES / csv), ["direct", "indirect", "rho_A"])
        write(f"tab_app_{tag}_est.tex", head_e, e)
        write(f"tab_app_{tag}_eq.tex", head_q, q)


def constants():
    """Design constants and application facts quoted in the text."""
    lines = []
    A = pd.read_csv(RES / "simA_size_power.csv")
    for func, cs in [("phi11_mean", "Phi"), ("direct", "Dir"), ("indirect", "Ind")]:
        lines.append(f"\\newcommand{{\\delta{cs}}}{{{A[A.functional == func].delta.iloc[0]:.4f}}}")
    F = pd.read_csv(RES / "app_data_facts.csv", index_col=0)["value"]
    lines.append(f"\\newcommand{{\\meanRain}}{{{F['mean_rain_mm']:.0f}}}")
    (TAB / "sim_constants.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("constants written")


if __name__ == "__main__":
    table_sim_summary("simA_size_power.csv", "tab_simA_summary.tex", [120, 240, 480])
    table_sim_summary("simC_both_misspecified.csv", "tab_simC_summary.tex", [120, 240, 480], nodeind=False, part="C")
    table_sim_summary("simD_application_size.csv", "tab_simD_summary.tex", [240, 480], beyond=1.5, nodeind=False, part="D")
    table_simCD()
    for func, lab in [("phi11_mean", "phi11"), ("direct", "direct"), ("indirect", "indirect")]:
        table_simA_full(func, lab)
    table_simB()
    if (RES / "simE_relative_margin.csv").exists():
        table_simE()
    table_app_similarity(); table_app_estimates(); table_app_eq12(); table_app_bench(); table_app_eq_rmse()
    table_app_cce(); table_app_jtest(); table_app_raw(); table_app_robust()
    constants()
