"""Candidate spatial weight matrices.

All matrices are row-standardised with a zero diagonal. The data-driven ones
(training correlation, graphical lasso) are built from the first `n_train`
months only, which in the application is the first rolling window, so they
never use a month that is later scored.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from sklearn.covariance import GraphicalLassoCV
from sklearn.preprocessing import StandardScaler


def row_norm(W):
    W = np.array(W, dtype=float, copy=True)
    np.fill_diagonal(W, 0.0)
    s = W.sum(axis=1, keepdims=True)
    s[s == 0] = 1.0
    return W / s


def build_knn(data_dir):
    """Geographic k-nearest-neighbour matrix (k=5) stored in data/knn_weights.csv."""
    A = pd.read_csv(f"{data_dir}/knn_weights.csv", index_col=0, float_precision="round_trip")
    return row_norm(A.to_numpy(dtype=float))


def haversine_km(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    dlat = la[:, None] - la[None, :]
    dlon = lo[:, None] - lo[None, :]
    hav = np.sin(dlat / 2) ** 2 + np.cos(la)[:, None] * np.cos(la)[None, :] * np.sin(dlon / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(hav, 0, 1)))


def build_invdist(data_dir, region_ids):
    cent = pd.read_csv(f"{data_dir}/regions.csv", float_precision="round_trip").set_index("region_id")
    lat = np.array([cent.loc[r, "lat"] for r in region_ids], float)
    lon = np.array([cent.loc[r, "lon"] for r in region_ids], float)
    dist = haversine_km(lat, lon)
    with np.errstate(divide="ignore"):
        W = np.where(dist > 0, 1.0 / dist, 0.0)
    return row_norm(W)


def build_corr(Y_train):
    """Positive sample correlations, row-standardised."""
    corr = np.corrcoef(Y_train.T)
    return row_norm(np.clip(corr, 0, None))


def build_glasso(Y_train, seed=0):
    """Positive partial correlations from a cross-validated graphical lasso."""
    Xs = StandardScaler().fit_transform(Y_train)
    model = GraphicalLassoCV().fit(Xs)
    P = model.precision_
    d = np.sqrt(np.diag(P))
    pcorr = -P / np.outer(d, d)
    return row_norm(np.clip(pcorr, 0, None))


def build_all(Y_train, data_dir, region_ids):
    return {
        "geographic KNN": build_knn(data_dir),
        "inverse-distance": build_invdist(data_dir, region_ids),
        "training-correlation": build_corr(Y_train),
        "graphical-lasso": build_glasso(Y_train),
    }
