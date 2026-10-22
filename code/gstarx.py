"""GSTARX(1;1) estimation, quantities of interest, and influence sequences.

Model, per location i (coefficients differ across locations):
    Y_{t,i} = c_i + phi10_i Y_{t-1,i} + phi11_i (W Y_{t-1})_i + gamma_i' X_{t-1} + u_{t,i}

Exogenous variables enter with one lag so that one-step forecasts only use
information available at t-1. All estimators are per-location OLS (identical to
the stacked block-diagonal OLS of the GSTAR literature).

Everything downstream treats OLS as estimating the pseudo-true (best linear
predictor) coefficient beta*(W), so nothing here assumes W is the true matrix.
"""
from __future__ import annotations
import numpy as np

# ----------------------------------------------------------------- design/OLS


def design(Y, W, X=None):
    """Regressors and targets.

    Y : (n, N) series (already seasonally adjusted); X : (n, q) or None.
    Returns Z (n-1, N, k) with z_{t,i} = (1, Y_{t-1,i}, (W Y_{t-1})_i, X_{t-1})
    and target (n-1, N) = Y[1:].
    """
    n, N = Y.shape
    F = Y @ W.T                                  # (n, N) spatial lag
    cols = [np.ones((n - 1, N)), Y[:-1], F[:-1]]
    if X is not None:
        for j in range(X.shape[1]):
            cols.append(np.repeat(X[:-1, j][:, None], N, axis=1))
    Z = np.stack(cols, axis=2)                   # (n-1, N, k)
    return Z, Y[1:]


def ols(Z, target):
    """Per-location OLS. Returns dict with beta (N,k), resid (n,N), Q (N,k,k),
    scores g (n,N,k) = z * u, and n."""
    Q = np.einsum("tik,til->ikl", Z, Z) / Z.shape[0]
    b = np.einsum("tik,ti->ik", Z, target) / Z.shape[0]
    beta = np.linalg.solve(Q, b[..., None])[..., 0]
    fitted = np.einsum("tik,ik->ti", Z, beta)
    resid = target - fitted
    g = Z * resid[..., None]
    return dict(beta=beta, resid=resid, Q=Q, g=g, n=Z.shape[0], fitted=fitted)


def fit(Y, W, X=None):
    Z, tgt = design(Y, W, X)
    out = ols(Z, tgt)
    out["Z"] = Z
    out["W"] = W
    return out

# ------------------------------------------------------- quantities of interest
# Column layout of beta: 0 = intercept, 1 = phi10, 2 = phi11, 3.. = gamma.


def phi11_mean(beta):
    return float(beta[:, 2].mean())


def phi11_mean_grad(beta):
    G = np.zeros_like(beta)
    G[:, 2] = 1.0 / beta.shape[0]
    return G


def transition(beta, W):
    return np.diag(beta[:, 1]) + np.diag(beta[:, 2]) @ W


def spectral_radius(A):
    return float(np.max(np.abs(np.linalg.eigvals(A))))


def impacts(beta, W, horizon=None):
    """Average direct, indirect and total effects of a one-time unit shock,
    cumulated over all horizons (horizon=None) or up to `horizon` months."""
    A = transition(beta, W)
    N = A.shape[0]
    I = np.eye(N)
    if horizon is None:
        M = np.linalg.solve(I - A, I)
    else:
        M = I.copy()
        P = I.copy()
        for _ in range(horizon):
            P = P @ A
            M = M + P
    tr = np.trace(M)
    tot = M.sum()
    return dict(direct=tr / N, indirect=(tot - tr) / N, total=tot / N,
                rho=spectral_radius(A), M=M)


def impact_grad(beta, W, which="direct"):
    """Gradient of the long-run average direct or indirect effect with respect
    to beta (N, k). Uses dM = M dA M with dA = E_ii dphi10_i and E_ii W dphi11_i.
      d tr(M)/d phi10_i = (M^2)_ii,            d tr(M)/d phi11_i = (W M^2)_ii
      d (1'M1)/d phi10_i = (1'M)_i (M1)_i,     d (1'M1)/d phi11_i = (1'M)_i (W M 1)_i
    """
    A = transition(beta, W)
    N = A.shape[0]
    M = np.linalg.solve(np.eye(N) - A, np.eye(N))
    M2 = M @ M
    dtr_10 = np.diag(M2)
    dtr_11 = np.diag(W @ M2)
    rowsum = M.sum(axis=0)          # (1'M)_i
    colsum = M.sum(axis=1)          # (M1)_i
    Wcol = W @ colsum               # (W M 1)_i
    dtot_10 = rowsum * colsum
    dtot_11 = rowsum * Wcol
    G = np.zeros_like(beta)
    if which == "direct":
        G[:, 1], G[:, 2] = dtr_10 / N, dtr_11 / N
    elif which == "indirect":
        G[:, 1], G[:, 2] = (dtot_10 - dtr_10) / N, (dtot_11 - dtr_11) / N
    elif which == "total":
        G[:, 1], G[:, 2] = dtot_10 / N, dtot_11 / N
    else:
        raise ValueError(which)
    return G


FUNCTIONALS = {
    "phi11_mean": (lambda f: phi11_mean(f["beta"]),
                   lambda f: phi11_mean_grad(f["beta"])),
    "direct": (lambda f: impacts(f["beta"], f["W"])["direct"],
               lambda f: impact_grad(f["beta"], f["W"], "direct")),
    "indirect": (lambda f: impacts(f["beta"], f["W"])["indirect"],
                 lambda f: impact_grad(f["beta"], f["W"], "indirect")),
}

# ------------------------------------------------------------ influence values


def influence_single(fitW, grad):
    """psi_t for tau(W) = h(beta(W)): sum_i grad_i' Q_i^{-1} g_{t,i}.  Returns (n,)."""
    Qinv_grad = np.linalg.solve(fitW["Q"], grad[..., None])[..., 0]   # (N, k)
    return np.einsum("tik,ik->t", fitW["g"], Qinv_grad)


def theta_and_influence(fitA, fitB, name):
    """theta_hat = h(beta_A) - h(beta_B) and its estimated influence sequence."""
    h, gradf = FUNCTIONALS[name]
    th = h(fitA) - h(fitB)
    psi = influence_single(fitA, gradf(fitA)) - influence_single(fitB, gradf(fitB))
    return float(th), psi

# ------------------------------------------------------- rolling-window forecast


def rolling_forecast(Yraw, month, W, X=None, m=240, score_idx=None, clim_fn=None,
                     anomaly_fn=None):
    """One-step forecasts with a fixed estimation window of m months.

    For each target month t = m, ..., n-1: climatology and coefficients are
    computed from months t-m .. t-1 only; the forecast of Y_t is the anomaly
    forecast plus the climatology of calendar month(t). Returns per-month loss
    (averaged over scored locations), per-node squared errors, and forecasts, all
    in the original units.
    """
    from data_io import climatology, anomalies
    n, N = Yraw.shape
    if score_idx is None:
        score_idx = np.arange(N)
    targets = np.arange(m, n)
    sq = np.empty((targets.size, score_idx.size))
    pred_all = np.empty((targets.size, N))
    for j, t in enumerate(targets):
        rows = np.arange(t - m, t)
        clim = climatology(Yraw, month, rows)
        Ya = anomalies(Yraw[rows], month[rows], clim)        # (m, N)
        Xw = X[rows] if X is not None else None
        f = fit(Ya, W, Xw)
        # regressor for forecasting Y_t uses Y_{t-1} (last window row) and X_{t-1}
        y_prev = Ya[-1]
        z = [np.ones(N), y_prev, W @ y_prev]
        if X is not None:
            for q in range(X.shape[1]):
                z.append(np.repeat(X[t - 1, q], N))
        z = np.stack(z, axis=1)                                # (N, k)
        pred_anom = np.einsum("ik,ik->i", z, f["beta"])
        pred = pred_anom + clim[month[t] - 1]
        pred_all[j] = pred
        sq[j] = (Yraw[t, score_idx] - pred[score_idx]) ** 2
    loss = sq.mean(axis=1)
    return dict(loss=loss, sq=sq, pred=pred_all, targets=targets,
                rmse=float(np.sqrt(sq.mean())))
