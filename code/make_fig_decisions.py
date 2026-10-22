"""Figure: boundary margins c_min (coefficients/effects) and RMSE intervals (forecasts)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
OUT = HERE.parent / "figures"; OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#b9b8b2"})
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
BLUE, ORANGE, AQUA, YEL = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
ab = {"geographic KNN": "KNN", "inverse-distance": "Inv. dist.", "training-correlation": "Corr.", "graphical-lasso": "G-lasso"}
Q = pd.read_csv(RES / "app_equivalence_Q12.csv")
R = pd.read_csv(RES / "app_equivalence_rmse.csv")
pairs = [f"{ab[a]} vs {ab[b]}" for a, b in zip(Q[Q.functional == "direct"].W_A, Q[Q.functional == "direct"].W_B)]
series = [("direct effect (total)", Q[Q.functional == "direct"].c_min_ref.values, BLUE, "o"),
          ("direct effect (W-dependent part)", Q[Q.functional == "direct"].c_min_dyn.values, AQUA, "s"),
          (r"mean spatial-lag coefficient $\bar\phi_{11}$", Q[Q.functional == "phi11_mean"].c_min_ref.values, YEL, "D"),
          ("indirect effect", Q[Q.functional == "indirect"].c_min_ref.values, ORANGE, "^")]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1.25, 1], wspace=0.55))
y = np.arange(len(pairs))[::-1]
off = np.linspace(-0.27, 0.27, len(series))
for (lab, v, c, m), o in zip(series, off):
    a1.scatter(100 * v, y + o, color=c, marker=m, s=22, label=lab, zorder=3, edgecolor="white", linewidth=0.5)
a1.axvline(10, color=INK, lw=1, ls=(0, (4, 3)))
a1.text(10.8, -0.75, "10% margin", fontsize=7.5, color=INK2)
a1.set_xscale("log"); a1.set_xlim(0.5, 300); a1.set_ylim(-1.0, len(pairs) - 0.5)
a1.set_xticks([1, 3, 10, 30, 100, 300]); a1.set_xticklabels(["1", "3", "10", "30", "100", "300"])
a1.minorticks_off()
a1.set_yticks(y); a1.set_yticklabels(pairs, fontsize=8)
a1.set_xlabel(r"smallest margin $c_{\min}$ (% of the reference value under KNN)")
a1.grid(axis="x", color=GRID); a1.set_axisbelow(True)
for s in ("top", "right"): a1.spines[s].set_visible(False)
a1.set_title("(a) coefficients and effects", loc="left", fontsize=9.5)
a1.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, handletextpad=0.3, columnspacing=1)
for k, (mdl, c, lab) in enumerate([("GSTAR", "#8a8984", "without climate indices"), ("GSTARX", BLUE, "with Niño 3.4 and DMI")]):
    r = R[R.model == mdl].reset_index(drop=True)
    o = -0.15 if k == 0 else 0.15
    for i in range(len(r)):
        yy = y[i] + o
        a2.plot([r.ci90_lo[i], r.ci90_hi[i]], [yy, yy], color=c, lw=2.2, solid_capstyle="round")
        a2.plot(r.theta_rmse_mm[i], yy, "o", color=c, ms=4, mec="white", mew=0.6)
    a2.plot([], [], color=c, lw=2.2, label=lab)
a2.axvspan(-1, 1, color="#e3eefb", zorder=0)
for v in (-1, 1): a2.axvline(v, color=BLUE, lw=0.9, ls=(0, (4, 3)))
a2.axvline(0, color=GRID, lw=1)
a2.set_yticks(y); a2.set_yticklabels(pairs, fontsize=8)
a2.set_xlim(-2.1, 2.1); a2.set_xlabel("RMSE difference (mm per month), 90% interval")
for s in ("top", "right"): a2.spines[s].set_visible(False)
a2.set_title(r"(b) one-step forecasts, margin $\pm1$ mm", loc="left", fontsize=9.5)
a2.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, frameon=False)
fig.savefig(OUT / "fig_decisions.pdf", bbox_inches="tight"); fig.savefig(OUT / "fig_decisions.png", dpi=150, bbox_inches="tight")
