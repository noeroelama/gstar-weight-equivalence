# gstar-weight-equivalence

Code and data for the paper

> A. Nurwahid. *A block bootstrap equivalence test for the choice of spatial weight matrix in
> space-time autoregressive models.* Submitted to *Computational Statistics & Data Analysis*.

The paper proposes a test of whether two candidate spatial weight matrices give the same value
of a reported quantity in a GSTARX(1;1) model, within a margin fixed in advance. The quantities
are the average spatial-lag coefficient, the long-run direct and indirect effects, and the
one-step forecast RMSE. The test is a two one-sided tests (TOST) procedure; its standard error
comes from a moving-block bootstrap of one scalar influence sequence per month, so the model is
not refitted inside the bootstrap. The application uses monthly CHIRPS rainfall for 40 regions
of Indonesia and its neighbours, January 1985 to April 2025.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | Python modules and scripts (see below) |
| `data/` | rainfall panel, climate indices, region table, nearest-neighbour matrix; sources in `data/README.md` |
| `results/` | the output files used in the paper (reference copies) |
| `figures/` | written by the figure scripts |

Core modules: `data_io.py` (data and monthly anomalies), `weights.py` (the four candidate
matrices), `gstarx.py` (per-location OLS, long-run effects and their gradients, rolling
forecasts), `inference.py` (influence sequences, vector moving-block bootstrap with the
Politis-White block length, Bartlett HAC, TOST with fixed and relative margins), `dgp.py`
(simulation designs and closed-form pseudo-true coefficients).

## Installation

Python 3.10 or later.

```
pip install -r requirements.txt
```

The results in the paper were computed with Python 3.10, NumPy 2.2, SciPy 1.15, pandas 2.3,
scikit-learn 1.7 and matplotlib 3.10.
With these versions the application scripts reproduce the files in `results/` up to rounding
error (differences below 1e-13). With newer versions (tested with NumPy 2.5, SciPy 1.18, pandas 3.0 and
scikit-learn 1.9) the graphical-lasso matrix changes slightly, some application numbers change
in the second or third decimal, and no equivalence decision changes.

## Reproducing the paper

Run the scripts from the `code/` folder. Each script writes to `../results/` and overwrites the
reference copies there.

| Script | Paper | Output |
|---|---|---|
| `python test_core.py` | checks of gradients, pseudo-true coefficients and standard errors | screen |
| `python app_run.py` | Section 6: Tables 3, 4 and 5 | `app_*.csv`, `app_summary.md` |
| `python app_extra.py` | Section 6: Tables 6 and 7; robustness tables in the Supplementary Material | `extra_*.csv`, `extra_summary.md` |
| `python app_revision.py` | Section 6: model order and regional effects; Table S14 | `rev_*.csv`, `rev_W.npz` |
| `python sim_study.py --exp AB` | Section 5: Parts A and B; Table 1 | `simA_size_power.csv`, `simB_coverage.csv` |
| `python sim_study.py --exp CD` | Section 5: Parts C and D; Table 2 | `simC_both_misspecified.csv`, `simD_application_size.csv` |
| `python sim_study.py --exp E` | Section 5: Part E (relative margin) | `simE_relative_margin.csv` |
| `python sim_boundary_reps.py --parts AC` and `--parts D` | per-replication records at the boundary (median standard errors in Tables 1 and 2) | `sim_boundary_reps_*.csv` |
| `python sim_lrv_check.py` | Section 5: long-run variance diagnostic | `sim_lrv_check.txt` |
| `python make_tables.py` | all tables as LaTeX | `results/tables/*.tex` |
| `python make_fig_sim.py` | Figure 1 | `figures/fig_simA.pdf` |
| `python make_fig_decisions.py` | Figure 3 | `figures/fig_decisions.pdf` |
| `python make_fig_spatial.py --provinces P.geojson --natural-earth NE.geojson` | Figures 2 and 4 (needs the polygons, see below) | `figures/fig_map_W.pdf`, `figures/fig_spill.pdf` |

The simulations use R = 1000 replications and B = 399 bootstrap draws with fixed seeds. The
runs reported in the paper took about 80 minutes each for `--exp AB` and `--exp CD` and about
5 minutes for `--exp E`. Smaller runs are possible with `--R` and `--B`. The application scripts
use B = 1999 and each takes less than a minute; `test_core.py` takes about 30 seconds.

## Polygons for the maps

`make_fig_spatial.py` needs two GeoJSON files that are not included here: the province layer
of https://github.com/Alf-Anas/batas-administrasi-indonesia (we used a simplified version) and
the Natural Earth 1:50m admin-0 countries file, for example
https://github.com/nvkelso/natural-earth-vector/blob/master/geojson/ne_50m_admin_0_countries.geojson.
The polygons are used only for drawing; all computations use `data/regions.csv` and
`data/knn_weights.csv`.

## Citation

If you use this code or data, please cite the paper and this archive (see `CITATION.cff`).

## License

Code: MIT License (see `LICENSE`). Data in `data/`: see `data/README.md` (rainfall series under
CC BY 4.0, derived from CHIRPS v2.0, which is in the public domain).
