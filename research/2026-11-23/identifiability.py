"""Linearized identifiability and conditional interval-coverage experiments."""
import numpy as np

def _validate(x, sigma):
    x, sigma = np.asarray(x, float), np.asarray(sigma, float)
    if x.ndim != 2 or sigma.shape != (len(x),) or len(x) <= x.shape[1]:
        raise ValueError("invalid design shape")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("nonfinite or nonpositive input")
    return x / sigma[:, None]

def singular_diagnostics(x, sigma):
    """Singular spectrum of whitened design; small values flag degeneracy."""
    a = _validate(x, sigma)
    s = np.linalg.svd(a, compute_uv=False)
    return s, float(np.inf if s[-1] == 0 else s[0] / s[-1])

def linear_estimate(x, y, sigma):
    """Full-rank conditional WLS; known independent Gaussian errors."""
    a = _validate(x, sigma)
    y = np.asarray(y, float)
    if y.shape != (len(a),) or not np.all(np.isfinite(y)):
        raise ValueError("invalid observations")
    q, r = np.linalg.qr(a, mode="reduced")
    if np.linalg.matrix_rank(r) < r.shape[0]:
        raise ValueError("unidentifiable coefficients")
    beta = np.linalg.solve(r, q.T @ (y / sigma))
    invr = np.linalg.solve(r, np.eye(r.shape[0]))
    covariance = invr @ invr.T
    return beta, covariance

def gaussian_coverage(x, sigma, truth, repetitions=400, seed=0, noise_scale=1.0):
    """Monte Carlo marginal 95% coverage against formal conditional intervals."""
    a = _validate(x, sigma)
    truth = np.asarray(truth, float)
    if truth.shape != (a.shape[1],) or repetitions < 2 or noise_scale <= 0:
        raise ValueError("invalid simulation parameters")
    if not np.isfinite(noise_scale) or np.linalg.matrix_rank(a) < a.shape[1]:
        raise ValueError("invalid scale or rank")
    rng = np.random.default_rng(seed)
    estimates = []
    standard_errors = None
    for _ in range(repetitions):
        y = np.asarray(x) @ truth + rng.normal(size=len(a)) * sigma * noise_scale
        beta, cov = linear_estimate(x, y, sigma)
        estimates.append(beta)
        standard_errors = np.sqrt(np.diag(cov))
    z = (np.asarray(estimates) - truth) / standard_errors
    return np.mean(np.abs(z) <= 1.959963984540054, axis=0), np.std(z, axis=0, ddof=1)
