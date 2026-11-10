"""Independent numerical diagnostics for conditional spectral-template inference.

No abundance interpretation; assumes supplied design matrix and covariance.
"""
import numpy as np


def gls_fit(design, observations, covariance):
    """Stable GLS via Cholesky whitening and SVD least squares.

    Raises ValueError on non-SPD covariance or rank-deficient design.
    """
    x = np.asarray(design, float)
    y = np.asarray(observations, float)
    c = np.asarray(covariance, float)
    if x.ndim != 2 or y.shape != (x.shape[0],) or c.shape != (len(y), len(y)):
        raise ValueError('incompatible dimensions')
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or not np.all(np.isfinite(c)):
        raise ValueError('nonfinite input')
    if not np.allclose(c, c.T, rtol=1e-10, atol=1e-12):
        raise ValueError('covariance must be symmetric')
    try:
        l = np.linalg.cholesky(c)
    except np.linalg.LinAlgError as exc:
        raise ValueError('covariance must be positive definite') from exc
    a = np.linalg.solve(l, x)
    b = np.linalg.solve(l, y)
    u, s, vt = np.linalg.svd(a, full_matrices=False)
    if not len(s) or s[-1] <= s[0] * max(a.shape) * np.finfo(float).eps:
        raise ValueError('design is rank deficient')
    beta = vt.T @ ((u.T @ b) / s)
    cov_beta = (vt.T / s**2) @ vt
    residual = b - a @ beta
    return beta, cov_beta, float(residual @ residual), len(y) - x.shape[1]


def ar1_covariance(sigma, rho):
    """Heteroscedastic AR(1) correlation, ordered by spectral bin index."""
    s = np.asarray(sigma, float)
    if s.ndim != 1 or len(s) < 2 or np.any(s <= 0) or not np.all(np.isfinite(s)):
        raise ValueError('sigma must be positive finite vector')
    if not np.isfinite(rho) or abs(rho) >= 1:
        raise ValueError('require |rho| < 1')
    d = np.abs(np.subtract.outer(np.arange(len(s)), np.arange(len(s))))
    return np.outer(s, s) * rho**d


def profile_template(design_nuisance, template, y, covariance):
    """GLS amplitude with nuisance terms profiled, including its 1-sigma error."""
    t = np.asarray(template, float)
    n = np.asarray(design_nuisance, float)
    if n.ndim != 2 or t.shape != (n.shape[0],):
        raise ValueError('bad template shape')
    x = np.column_stack((n, t))
    beta, cov, chi2, dof = gls_fit(x, y, covariance)
    return float(beta[-1]), float(np.sqrt(cov[-1, -1])), chi2, dof
