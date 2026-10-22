"""Figure 1 of Section 5: share of replications declaring equivalence against the true
difference |theta|/delta, Part A, one panel per quantity, one line per sample length
(vector MBB). Reads ../results/simA_size_power.csv; writes ../figures/fig_simA.pdf.
Colours: the validated categorical palette used for the graphical abstract.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
OUT = HERE.parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})

A = pd.read_csv(RES / "simA_size_power.csv")
panels = [("phi11_mean", r"average spatial-lag coefficient $\bar\phi_{11}$"),
          ("direct", "long-run average direct effect"),
          ("indirect", "long-run average indirect effect")]
fig, axes = plt.subplots(1, 3, figsize=(16 / 2.54, 5.6 / 2.54), sharey=True)
for ax, (func, title) in zip(axes, panels):
    D = A[A.functional == func]
    for T, col, ls in [(120, AQUA, (0, (1.2, 1.2))), (240, ORANGE, (0, (4, 1.5))), (480, BLUE, "-")]:
        g = D[D["T"] == T].sort_values("ratio")
        ax.plot(g.ratio, g.eq_mbb_pw, color=col, lw=1.5, ls=ls, marker="o", ms=2.6, label=f"T = {T}")
    ax.axvline(1.0, color=INK2, lw=0.6)
    ax.axhline(0.05, color=INK2, lw=0.6, ls=(0, (2, 2)))
    ax.set_title(title, fontsize=7.5, color=INK, loc="left", pad=3)
    ax.set_xlim(0, 2.05); ax.set_ylim(-0.02, 1.03)
    ax.set_xticks([0, 0.5, 1, 1.5, 2]); ax.set_yticks([0.05, 0.5, 1])
    ax.set_yticklabels(["0.05", "0.5", "1"])
    ax.grid(axis="y", color=GRID, lw=0.5)
    ax.tick_params(length=2, pad=1.5)
    ax.set_xlabel(r"true difference $|\theta|/\delta$", labelpad=1)
axes[0].set_ylabel("share declaring equivalence", labelpad=2)
axes[0].text(1.04, 0.62, "boundary\n" r"$|\theta|=\delta$", fontsize=6.5, color=INK2)
axes[2].legend(frameon=False, fontsize=6.8, loc="upper right", handlelength=2.6)
fig.tight_layout(pad=0.4, w_pad=0.8)
fig.savefig(OUT / "fig_simA.pdf")
fig.savefig(OUT / "fig_simA.png", dpi=200)
print("saved", OUT / "fig_simA.pdf")
