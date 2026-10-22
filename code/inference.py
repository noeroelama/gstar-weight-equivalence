"""Standard errors for a mean, and the TOST equivalence decision.

Three estimators of the long-run variance of a stationary sequence v_t:
  * mbb_se        : moving-block bootstrap (Kuensch 1989) on the observed sequence.
                    When v_t is the sum over locations of per-location terms, the
                    sequence already contains the whole cross-section, so this is
                    the vector MBB of the paper.
  * node_indep_se : the comparator that resamples each location's contribution
                    with its own independent block draws (breaks cross-sectional
                    dependence). Used only to show what goes wrong.
  * hac_se        : Bartlett-kernel (Newey-West) long-run variance.
Plus the Politis-White automatic block length and the TOST decision.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import norm

# ------------------------------------------------------------ block length


def politis_white(x):
    """Politis & White (2004) with the Patton-Politis-White (2009) correction.
    Returns the block length for the moving-block bootstrap (b_cb)."""
    x = np.asarray(x, float)
    n = x.size
    eps = x - x.mean()
    kn = max(5, int(np.ceil(np.log10(n))))
    mmax = int(np.ceil(np.sqrt(n))) + kn
    bmax = np.ceil(min(3.0 * np.sqrt(n), n / 3.0))
    acv = np.array([np.dot(eps[: n - k], eps[k:]) / n for k in range(mmax + 1)])
    rho = acv / acv[0]
    cv = 2.0 * np.sqrt(np.log10(n) / n)
    insig = np.abs(rho) < cv
    m_hat = None
    for m in range(1, mmax - kn + 1):
        if np.all(insig[m + 1: m + kn + 1]):
            m_hat = m
            break
    if m_hat is None:
        sig = np.where(np.abs(rho[1:]) >= cv)[0]
        m_hat = int(sig[-1] + 1) if sig.size else 1
    M = max(1, min(2 * m_hat, mmax))
    kk = np.arange(-M, M + 1)
    at = np.abs(kk / M)
    lam = np.where(at <= 0.5, 1.0, np.where(at <= 1.0, 2.0 * (1.0 - at), 0.0))
    acv_sym = acv[np.abs(kk)]
    g_hat = np.sum(lam * np.abs(kk) * acv_sym)
    G = np.sum(lam * acv_sym)
    d_cb = (4.0 / 3.0) * G ** 2
    b_cb = ((2.0 * g_hat ** 2) / d_cb) ** (1 / 3) * n ** (1 / 3)
    return float(max(1.0, min(b_cb, bmax)))

# ------------------------------------------------------------ bootstrap


def block_length(L):
    """Block length actually used: the automatic length rounded up, and at least 2, so
    that the bootstrap never degenerates to resampling single periods."""
    return int(max(2, np.ceil(L)))


def block_index(n, L, B, rng):
    """(B, n) indices: ceil(n/L) blocks of L consecutive positions, with replacement.
    L is passed through block_length (rounded up, at least 2)."""
    L = int(min(block_length(L), n))
    nb = int(np.ceil(n / L))
    starts = rng.integers(0, n - L + 1, size=(B, nb))
    idx = starts[:, :, None] + np.arange(L)[None, None, :]
    return idx.reshape(B, nb * L)[:, :n]


def mbb_means(v, L, B, rng):
    """Bootstrap means of a scalar sequence v (n,). Returns (B,)."""
    idx = block_index(v.size, L, B, rng)
    return v[idx].mean(axis=1)


def mbb_se(v, L, B, rng):
    """Standard error of the sample mean of v via the moving-block bootstrap."""
    return float(mbb_means(np.asarray(v, float), L, B, rng).std(ddof=1))


def node_indep_se(contrib, L, B, rng):
    """contrib: (n, N) per-location contributions with v_t = contrib[t].sum().
    Each location gets its own block draws, so the cross-section is broken."""
    n, N = contrib.shape
    means = np.zeros(B)
    for i in range(N):
        idx = block_index(n, L, B, rng)
        means += contrib[idx, i].mean(axis=1)
    return float(means.std(ddof=1))

# ------------------------------------------------------------ HAC


def nw_bandwidth(n):
    """Newey-West (1994) style rule of thumb: floor(4 (n/100)^(2/9))."""
    return int(np.floor(4.0 * (n / 100.0) ** (2.0 / 9.0)))


def hac_lrv(v, bandwidth=None):
    """Bartlett-kernel long-run variance of a scalar sequence."""
    v = np.asarray(v, float)
    n = v.size
    e = v - v.mean()
    if bandwidth is None:
        bandwidth = nw_bandwidth(n)
    s = np.dot(e, e) / n
    for k in range(1, bandwidth + 1):
        w = 1.0 - k / (bandwidth + 1.0)
        s += 2.0 * w * np.dot(e[:-k], e[k:]) / n
    return float(max(s, 0.0))


def hac_se(v, bandwidth=None):
    v = np.asarray(v, float)
    return float(np.sqrt(hac_lrv(v, bandwidth) / v.size))

# ------------------------------------------------------------ TOST


def tost(theta_hat, se, delta, alpha=0.05):
    """Two one-sided tests. Equivalent iff the (1-2alpha) interval is inside (-delta, delta).
    Also returns the smallest delta at which equivalence would be declared."""
    z = norm.ppf(1 - alpha)
    lo, hi = theta_hat - z * se, theta_hat + z * se
    t1 = (theta_hat + delta) / se        # H01: theta <= -delta, reject if t1 > z
    t2 = (theta_hat - delta) / se        # H02: theta >= +delta, reject if t2 < -z
    p1 = 1 - norm.cdf(t1)
    p2 = norm.cdf(t2)
    return dict(theta=theta_hat, se=se, ci_lo=lo, ci_hi=hi,
                equivalent=bool(lo > -delta and hi < delta),
                p_tost=float(max(p1, p2)),
                min_delta=float(max(abs(lo), abs(hi))))

# ------------------------------------------------------------ relative margin


def hac_cov2(a, b, bandwidth=None):
    """Bartlett-kernel long-run covariance matrix of (a_t, b_t), divided by n, so that it
    estimates the covariance matrix of the two sample means."""
    Zm = np.column_stack([a, b]).astype(float)
    n = Zm.shape[0]
    E = Zm - Zm.mean(axis=0)
    if bandwidth is None:
        bandwidth = nw_bandwidth(n)
    S = E.T @ E / n
    for k in range(1, bandwidth + 1):
        w = 1.0 - k / (bandwidth + 1.0)
        Gk = E[k:].T @ E[:-k] / n
        S += w * (Gk + Gk.T)
    return S / n


def tost_relative(theta_hat, psi, tau_ref, psi_ref, c=0.10, B=1999, rng=None,
                  alphas=(0.05, 0.025), c_max=5.0, V=None):
    """Equivalence with a margin proportional to a reference quantity:
        H0: |theta| >= c * tau_ref   vs   H1: |theta| < c * tau_ref,   tau_ref > 0.
    H0 is the union of {theta + c tau_ref <= 0} and {theta - c tau_ref >= 0}. Each
    one-sided hypothesis is a smooth function of the pseudo-true coefficients, with
    influence sequence psi +/- c psi_ref, so it is tested with its own standard error.
    The two sequences are bootstrapped jointly with the same block indices, which gives
    se(c, +/-) = sqrt(V11 +/- 2 c V12 + c^2 V22) for every c from one 2x2 bootstrap
    covariance. Returns the decision at c, the smallest c at which equivalence is
    declared (c_min) for each alpha, and the standard errors."""
    rng = np.random.default_rng() if rng is None else rng
    psi, psi_ref = np.asarray(psi, float), np.asarray(psi_ref, float)
    n = psi.size
    if V is None:                       # vector MBB, the same blocks for both sequences
        Lhat = max(politis_white(psi), politis_white(psi_ref))
        idx = block_index(n, Lhat, B, rng)
        m = np.column_stack([psi[idx].mean(axis=1), psi_ref[idx].mean(axis=1)])
        V = np.cov(m, rowvar=False, ddof=1)
    else:                               # e.g. hac_cov2(psi, psi_ref)
        Lhat = float("nan")
    se = lambda cc, sgn: float(np.sqrt(max(V[0, 0] + sgn * 2 * cc * V[0, 1] + cc ** 2 * V[1, 1], 0.0)))
    out = dict(theta=float(theta_hat), tau_ref=float(tau_ref), L_pw=float(Lhat),
               L_used=block_length(Lhat) if np.isfinite(Lhat) else None, se_theta=float(np.sqrt(V[0, 0])),
               se_ref=float(np.sqrt(V[1, 1])), cov=float(V[0, 1]))
    grid = np.linspace(0.0, c_max, int(c_max / 0.0005) + 1)
    lo = theta_hat + grid * tau_ref          # must be significantly > 0
    hi = theta_hat - grid * tau_ref          # must be significantly < 0
    se_p = np.sqrt(np.maximum(V[0, 0] + 2 * grid * V[0, 1] + grid ** 2 * V[1, 1], 0))
    se_m = np.sqrt(np.maximum(V[0, 0] - 2 * grid * V[0, 1] + grid ** 2 * V[1, 1], 0))
    for a in alphas:
        z = norm.ppf(1 - a)
        ok = (lo - z * se_p > 0) & (hi + z * se_m < 0)
        tag = f"{int(round(a * 1000)):03d}"
        out[f"c_min_{tag}"] = float(grid[np.argmax(ok)]) if ok.any() else float("nan")
        out[f"equivalent_{tag}"] = bool((theta_hat + c * tau_ref - z * se(c, +1) > 0) and
                                        (theta_hat - c * tau_ref + z * se(c, -1) < 0))
    return out

# ------------------------------------------------------------ RMSE difference


def rmse_diff_sequence(loss_a, loss_b):
    """theta_hat = sqrt(mean loss_a) - sqrt(mean loss_b) and the delta-method
    linearised sequence whose mean has the same limit distribution."""
    ma, mb = loss_a.mean(), loss_b.mean()
    theta = float(np.sqrt(ma) - np.sqrt(mb))
    d_tilde = (loss_a - ma) / (2 * np.sqrt(ma)) - (loss_b - mb) / (2 * np.sqrt(mb))
    return theta, d_tilde
