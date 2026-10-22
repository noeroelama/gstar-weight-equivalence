"""Load the CHIRPS panel and build seasonally adjusted series.

Data files in data/ (see data/README.md):
  rainfall_monthly.csv  date + one column per region, mm/month, 1985-01 to 2025-04 (T=484, N=40)
  climate_indices.csv   date, nino34, dmi
  regions.csv           region_id, country, name, kode_prov, lat, lon (same order as the columns)

Seasonal adjustment: monthly anomalies. For each location and calendar month,
subtract the mean over a stated set of months (the "climatology window"). In the
rolling-window forecast the climatology is computed from the same window that is
used for estimation, so no future information is used.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
RT = dict(float_precision="round_trip")   # read floats exactly as written


def load_panel(data_dir: Path | None = None):
    data_dir = Path(data_dir or DATA)
    reg = pd.read_csv(data_dir / "regions.csv", dtype={"kode_prov": str})
    region_ids = list(reg["region_id"])
    Ydf = pd.read_csv(data_dir / "rainfall_monthly.csv", **RT)
    Xdf = pd.read_csv(data_dir / "climate_indices.csv", **RT)
    assert list(Ydf["date"]) == list(Xdf["date"]), "dates of rainfall and climate files differ"
    Y = Ydf[region_ids].to_numpy(dtype=float)  # (T, N)
    X = Xdf[["nino34", "dmi"]].to_numpy(dtype=float)  # (T, 2)
    dates = Ydf["date"].to_numpy()
    countries = list(reg["country"])
    idn = np.array([i for i, c in enumerate(countries) if c == "IDN"])
    month = np.array([int(str(s)[5:7]) for s in dates])  # 1..12
    return dict(Y=Y, X=X, dates=dates, region_ids=region_ids,
                countries=countries, idn=idn, month=month)


def climatology(Y: np.ndarray, month: np.ndarray, rows: np.ndarray):
    """Per-location mean for each calendar month, using only the rows given.
    Returns (12, N)."""
    clim = np.zeros((12, Y.shape[1]))
    for mth in range(1, 13):
        sel = rows[month[rows] == mth]
        clim[mth - 1] = Y[sel].mean(axis=0)
    return clim


def anomalies(Y: np.ndarray, month: np.ndarray, clim: np.ndarray):
    """Y minus the climatology of its calendar month. Shapes: Y (n, N), month (n,)."""
    return Y - clim[month - 1]
