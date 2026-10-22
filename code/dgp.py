"""Data-generating process for the simulation study, with exact population
pseudo-true coefficients for any candidate W.

DGP (no exogenous variables):
    Y_t = Phi10 Y_{t-1} + Phi11 W0 Y_{t-1} + eps_t,   eps_t ~ N(0, Sigma) i.i.d.
with location-specific phi10_i, phi11_i and spatially correlated Sigma
(exponential in distance between random locations on the unit square).

Population quantities. With A0 = Phi10 + Phi11 W0 and Gamma0 = Var(Y_t) solving
Gamma0 = A0 Gamma0 A0' + Sigma, the regressor z_{t,i}(W) = (1, Y_{t-1,i}, (W Y_{t-1})_i)
has second-moment matrix built from Gamma0, and E[z_{t,i} Y_{t,i}] is built from
Gamma1 = Cov(Y_t, Y_{t-1}) = A0 Gamma0. Hence beta*_i(W) = Q_i(W)^{-1} E[z Y] is
available in closed form for every W, true or not.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import solve_discrete_lyapunov

from weights import row_norm


def knn_from_coords(coords, k):
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    W = np.zeros((len(coords), len(coords)))
    for i in range(len(coords)):
        W[i, np.argsort(d[i])[:k]] = 1.0
    return row_norm(W)


def band_from_coords(coords, mean_neighbours=4.0):
    """Distance-band graph: i and j are neighbours when their distance is below a radius
    chosen so that the average number of neighbours is about `mean_neighbours`; a location
    with no neighbour inside the band is linked to its nearest location. Row-normalised."""
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    N = len(coords)
    off = np.sort(d[np.isfinite(d)])
    radius = off[min(off.size - 1, int(round(mean_neighbours * N)) - 1)]
    W = (d <= radius).astype(float)
    for i in range(N):
        if W[i].sum() == 0:
            W[i, np.argmin(d[i])] = 1.0
    return row_norm(W)


def make_dgp(N=18, k=3, corr_range=0.6, phi10=(0.10, 0.30), phi11=(0.40, 0.60),
             sigma=1.0, seed=0, graph="knn", q=0, innov="gauss", df=5):
    """Random locations, true matrix W0 (KNN with k neighbours, or a distance band),
    heterogeneous coefficients, spatially correlated innovations. Optionally q exogenous
    variables X_t = Pi X_{t-1} + xi_t entering as gamma_i' X_{t-1}, and Student-t
    innovations (scaled to the same covariance). Rescales coefficients if the system
    is not stable."""
    rng = np.random.default_rng(seed)
    coords = rng.random((N, 2))
    W0 = knn_from_coords(coords, k) if graph == "knn" else band_from_coords(coords)
    p10 = rng.uniform(*phi10, size=N)
    p11 = rng.uniform(*phi11, size=N)
    A0 = np.diag(p10) + np.diag(p11) @ W0
    rho = np.max(np.abs(np.linalg.eigvals(A0)))
    if rho >= 0.95:
        scale = 0.9 / rho
        p10, p11 = p10 * scale, p11 * scale
        A0 = np.diag(p10) + np.diag(p11) @ W0
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
    Sigma = sigma ** 2 * np.exp(-d / corr_range) if corr_range > 0 else sigma ** 2 * np.eye(N)
    Gam = rng.uniform(-0.3, 0.3, size=(N, q)) if q > 0 else np.zeros((N, 0))
    Pi = np.diag(np.linspace(0.8, 0.6, q)) if q > 0 else np.zeros((0, 0))
    # stacked state S_t = (Y_t, X_t): S_t = A_S S_{t-1} + (eps_t, xi_t)
    A_S = np.block([[A0, Gam], [np.zeros((q, N)), Pi]])
    Sig_S = np.block([[Sigma, np.zeros((N, q))], [np.zeros((q, N)), np.eye(q)]])
    G0 = solve_discrete_lyapunov(A_S, Sig_S)
    return dict(N=N, q=q, coords=coords, W0=W0, phi10=p10, phi11=p11, A0=A0, Gam=Gam, Pi=Pi,
                Sigma=Sigma, L=np.linalg.cholesky(Sigma), innov=innov, df=df,
                A_S=A_S, G0S=G0, G1S=A_S @ G0,
                Gamma0=G0[:N, :N], Gamma1=(A_S @ G0)[:N, :N],
                rho=np.max(np.abs(np.linalg.eigvals(A0))))


def simulate(dgp, T, rng, burn=200):
    """Returns Y (T, N), or (Y, X) when the DGP has exogenous variables."""
    N, q = dgp["N"], dgp.get("q", 0)
    if dgp.get("innov", "gauss") == "t":
        df = dgp["df"]
        Z = rng.standard_t(df, size=(T + burn, N)) * np.sqrt((df - 2) / df)   # unit variance
    else:
        Z = rng.standard_normal((T + burn, N))
    E = Z @ dgp["L"].T
    if q == 0:
        Y = np.zeros((T + burn, N))
        A0 = dgp["A0"]
        for t in range(1, T + burn):
            Y[t] = A0 @ Y[t - 1] + E[t]
        return Y[burn:]
    Xi = rng.standard_normal((T + burn, q))
    S = np.zeros((T + burn, N + q))
    A_S = dgp["A_S"]
    for t in range(1, T + burn):
        S[t] = A_S @ S[t - 1] + np.concatenate([E[t], Xi[t]])
    return S[burn:, :N], S[burn:, N:]


def pseudo_true(dgp, W):
    """beta*_i(W) = (c, phi10, phi11, gamma) for each location, in closed form."""
    if dgp.get("q", 0) > 0:
        return _pseudo_true_x(dgp, W)
    G0, G1 = dgp["Gamma0"], dgp["Gamma1"]
    N = dgp["N"]
    beta = np.zeros((N, 3))
    # second moments of (Y_{t-1,i}, (W Y_{t-1})_i) and cross moments with Y_{t,i}
    WG0 = W @ G0                      # Cov(W Y, Y)
    WG0W = W @ G0 @ W.T               # Var(W Y)
    G1W = G1 @ W.T                    # Cov(Y_t, W Y_{t-1})
    for i in range(N):
        Q = np.array([[1.0, 0.0, 0.0],
                      [0.0, G0[i, i], WG0[i, i]],
                      [0.0, WG0[i, i], WG0W[i, i]]])
        b = np.array([0.0, G1[i, i], G1W[i, i]])
        beta[i] = np.linalg.solve(Q, b)
    return beta


def _pseudo_true_x(dgp, W):
    """Closed form with exogenous variables. z_{t,i} = (1, R_i S_{t-1}) with
    R_i = [e_i' ; (W_i, 0) ; (0, I_q)]; the process has mean zero, so the intercept is 0,
    Q_i = blockdiag(1, R_i Var(S) R_i') and E[z Y_{t,i}] = (0, R_i Cov(S_{t-1}, Y_{t,i}))."""
    N, q = dgp["N"], dgp["q"]
    G0, G1 = dgp["G0S"], dgp["G1S"]            # Var(S_t), E[S_t S_{t-1}']
    beta = np.zeros((N, 3 + q))
    for i in range(N):
        R = np.zeros((2 + q, N + q))
        R[0, i] = 1.0
        R[1, :N] = W[i]
        R[2:, N:] = np.eye(q)
        Qs = R @ G0 @ R.T
        b = R @ G1[i, :]                       # Cov(S_{t-1}, Y_{t,i})
        beta[i, 1:] = np.linalg.solve(Qs, b)
    return beta


def blend(WA, Wwrong, lam):
    """Candidate W_B(lam) moving from W_A (lam=0) towards a different graph."""
    return row_norm((1 - lam) * WA + lam * Wwrong)
