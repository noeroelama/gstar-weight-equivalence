"""Spatial figures: study area with the candidate matrices (Figure 2) and regional spillovers (Figure 4).

The region polygons are not distributed with this repository. Provide two GeoJSON files:
  --provinces      Indonesian provinces (we used the province layer of
                   https://github.com/Alf-Anas/batas-administrasi-indonesia, simplified)
  --natural-earth  Natural Earth 1:50m admin-0 countries
                   (ne_50m_admin_0_countries.geojson, https://www.naturalearthdata.com)
Usage: python make_fig_spatial.py --provinces PROV.geojson --natural-earth NE50.geojson
Reads ../data/regions.csv, ../results/rev_W.npz, ../results/rev_region_effects.csv;
writes ../figures/fig_map_W.pdf and ../figures/fig_spill.pdf.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--provinces", required=True, type=Path)
ap.add_argument("--natural-earth", required=True, type=Path)
args = ap.parse_args()

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"
DATA = HERE.parent / "data"
OUT = HERE.parent / "figures"; OUT.mkdir(parents=True, exist_ok=True)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"
BLUE, ORANGE = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#b9b8b2", "axes.titlesize": 9.5})

cen = pd.read_csv(DATA / "regions.csv", float_precision="round_trip").set_index("region_id")
Wz = np.load(RES / "rev_W.npz", allow_pickle=True)
rid = [str(r) for r in Wz["region_ids"]]
lon = np.array([cen.loc[r, "lon"] for r in rid]); lat = np.array([cen.loc[r, "lat"] for r in rid])
Wn = {"geographic KNN": Wz["geographic_KNN"], "graphical lasso": Wz["graphical_lasso"],
      "inverse distance": Wz["inverse_distance"], "training correlation": Wz["training_correlation"]}

# --- polygons, each assigned to the region whose centroid is nearest the polygon's centroid
def rings(g):
    return [np.array(g["coordinates"][0])[:, :2]] if g["type"] == "Polygon" else [np.array(p[0])[:, :2] for p in g["coordinates"]]
prov = json.load(open(args.provinces, encoding="utf-8"))
ne = json.load(open(args.natural_earth, encoding="utf-8"))
polys = {r: [] for r in rid}; context = []
def nearest(x, y, pool):
    d = (lon[pool] - x) ** 2 + (lat[pool] - y) ** 2
    return pool[int(np.argmin(d))]
idn_idx = np.array([i for i, r in enumerate(rid) if r.startswith("ID-")])
for f in prov["features"]:
    rs = rings(f["geometry"])
    big = max(rs, key=len)
    i = nearest(big[:, 0].mean(), big[:, 1].mean(), idn_idx)
    polys[rid[i]] += rs
for f in ne["features"]:
    nm = f["properties"].get("ADMIN")
    for r in rings(f["geometry"]):
        cx, cy = r[:, 0].mean(), r[:, 1].mean()
        if not (88 < cx < 160 and -16 < cy < 12) or nm == "Indonesia":
            continue
        if nm == "Papua New Guinea":
            polys["PNG_ALL"].append(r)
        elif nm == "Brunei" or (nm == "Malaysia" and cx > 108.5):
            polys["BORNEO_NONIDN"].append(r)
        else:
            context.append(r)
assert all(len(v) for v in polys.values()), [k for k, v in polys.items() if not v]

def base(ax, values=None, cmap=None, norm=None, fill="#e3eefb"):
    for r in context:
        cfc = "#f1f0ec" if values is None else "white"
        ax.add_patch(Polygon(r, closed=True, fc=cfc, ec="#c9c8c2", lw=0.3, zorder=0))
    for k, rs in polys.items():
        fc = fill if values is None else cmap(norm(values[k]))
        for r in rs:
            ax.add_patch(Polygon(r, closed=True, fc=fc, ec="white", lw=0.35, zorder=0))
    ax.set_rasterization_zorder(0.5)
    ax.set_xlim(94, 156); ax.set_ylim(-11.5, 8.5); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

def edges(ax, W, color, wmin=0.0):
    N = len(rid)
    for i in range(N):
        for j in range(i + 1, N):
            w = max(W[i, j], W[j, i])
            if w > wmin:
                ax.plot([lon[i], lon[j]], [lat[i], lat[j]], color=color, lw=0.4 + 2.2 * w, alpha=0.75, zorder=2)
    ax.scatter(lon, lat, s=7, color=INK, zorder=3)

# ------------------------------------------------------------- Figure: study area and W
order = np.argsort(lon)
fig = plt.figure(figsize=(7.2, 4.3))
gs = fig.add_gridspec(2, 4, height_ratios=[0.62, 1], hspace=0.12, wspace=0.12)
for c, (nm, col) in enumerate([("geographic KNN", BLUE), ("graphical lasso", BLUE)]):
    ax = fig.add_subplot(gs[0, 2 * c:2 * c + 2]); base(ax); edges(ax, Wn[nm], col)
    ax.set_title(f"({'ab'[c]}) {nm}: {(Wn[nm] > 0).sum(1).mean():.2g} links per row", loc="left")
vmax = max(W.max() for W in Wn.values())
for c, nm in enumerate(Wn):
    ax = fig.add_subplot(gs[1, c])
    im = ax.imshow(Wn[nm][np.ix_(order, order)], cmap=LinearSegmentedColormap.from_list("b", ["#ffffff", "#9cc3f0", BLUE, "#123d70"]),
                   vmin=0, vmax=0.5, interpolation="nearest")
    ax.set_title(f"({'cdef'[c]}) {nm}", loc="left", fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([])
    if c == 0:
        ax.set_ylabel("row (region, west to east)", fontsize=8)
    ax.set_xlabel("column", fontsize=8)
cax = fig.add_axes([0.92, 0.11, 0.012, 0.42]); cb = fig.colorbar(im, cax=cax); cb.set_label("weight", fontsize=8)
cb.ax.tick_params(labelsize=7)
fig.savefig(OUT / "fig_map_W.pdf", bbox_inches="tight", dpi=300); fig.savefig(OUT / "fig_map_W.png", dpi=150, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------- Figure: regional spillover received
eff = pd.read_csv(RES / "rev_region_effects.csv")
P = eff.pivot(index="region_id", columns="W", values="spill_received")
a, b = P["geographic KNN"], P["training-correlation"]
cm = LinearSegmentedColormap.from_list("div", [ORANGE, "#f6c2a8", "#efeeea", "#a9c9f2", BLUE])
lim = float(np.ceil(max(abs(P[["geographic KNN", "training-correlation"]].values).max(), 0.9) * 10) / 10)
norm = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
d = b - a
dl = float(np.ceil(abs(d).max() * 10) / 10)
dnorm = TwoSlopeNorm(vmin=-dl, vcenter=0, vmax=dl)
fig, axs = plt.subplots(3, 1, figsize=(6.4, 7.0))
for ax, vals, t, nr, cl in [(axs[0], a, "(a) geographic KNN (sparse): spillover received", norm, "spillover received"),
                        (axs[1], b, "(b) training correlation (dense): spillover received", norm, "spillover received"),
                        (axs[2], d, "(c) difference, (b) minus (a)", dnorm, "difference")]:
    base(ax, vals.to_dict(), cm, nr)
    ax.set_title(t, loc="left")
    sm = plt.cm.ScalarMappable(norm=nr, cmap=cm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01); cb.ax.tick_params(labelsize=7); cb.set_label(cl, fontsize=7.5)
fig.tight_layout()
fig.savefig(OUT / "fig_spill.pdf", bbox_inches="tight", dpi=300); fig.savefig(OUT / "fig_spill.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("mean received", a.mean().round(3), b.mean().round(3), "corr", np.corrcoef(a, b)[0, 1].round(2))
print("largest differences:\n", d.sort_values(ascending=False).round(2).head(8).to_string())
print("smallest:\n", d.sort_values().round(2).head(4).to_string())
