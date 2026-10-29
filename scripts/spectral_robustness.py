"""Independent, conditional linear diagnostics for the pinned WASP-39 b spectrum.

Gaussian features are descriptive basis functions, NOT molecular opacities.
The within-product AR(1) covariance is a hypothetical sensitivity model,
not a measured JWST covariance matrix. No abundance/detection inference.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def _gauss(w: np.ndarray, center: float, width: float) -> np.ndarray:
    return np.exp(-0.5 * ((w - center) / width) ** 2)


def templates(w: np.ndarray, width_scale: float = 1.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Match the public browser's reduced-order template definitions."""
    if not np.isfinite(width_scale) or width_scale <= 0:
        raise ValueError("width_scale must be finite and positive")
    g = lambda c, s: _gauss(w, c, s * width_scale)
    return (0.62 * g(1.15, .08) + g(1.40, .12) + .42 * g(1.90, .14)
            + .48 * g(2.70, .18), g(4.30, .16), g(4.05, .07))


def fit_spectrum(points: list[dict], *, rho: float = 0.0,
                 error_floor_ppm: float = 0.0, width_scale: float = 1.0) -> dict:
    """Full-rank GLS with per-product offsets, slope and three templates.

    AR(1) correlation is by increasing wavelength *index* within each product;
    it is not a physical correlation length in micrometres. Returned errors
    are conditional on the assumed covariance and fixed template family.
    """
    if not np.isfinite(rho) or abs(rho) >= 1:
        raise ValueError("rho must satisfy |rho| < 1")
    if not np.isfinite(error_floor_ppm) or error_floor_ppm < 0:
        raise ValueError("error floor must be finite and nonnegative")
    if len(points) < 12:
        raise ValueError("at least twelve spectral bins required")
    w = np.array([p["wavelength_um"] for p in points], dtype=float)
    y = np.array([p["transit_depth_ppm"] for p in points], dtype=float)
    s = np.array([p["uncertainty_ppm"] for p in points], dtype=float)
    names = np.array([str(p["instrument"]) for p in points])
    if not (np.all(np.isfinite(w)) and np.all(w > 0) and
            np.all(np.isfinite(y)) and np.all(np.isfinite(s)) and np.all(s > 0)):
        raise ValueError("wavelength, depth and uncertainty must be finite and physical")
    instruments = sorted(set(names.tolist()))
    h2o, co2, so2 = templates(w, width_scale)
    columns = [f"offset:{name}" for name in instruments] + ["log_slope", "h2o", "co2", "so2"]
    x = np.column_stack([*(names == name for name in instruments),
                         np.log(w / 2.2), h2o, co2, so2]).astype(float)
    if len(points) <= x.shape[1] + 1:
        raise ValueError("not enough degrees of freedom")
    sigma = np.hypot(s, error_floor_ppm)
    whitened_x, whitened_y = [], []
    for name in instruments:
        idx = np.flatnonzero(names == name)
        idx = idx[np.argsort(w[idx], kind="stable")]
        n = len(idx)
        # Cholesky whitening of D R D with R_ij = rho**abs(i-j).
        distance = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
        correlation = rho ** distance
        chol = np.linalg.cholesky(correlation)
        whitened_x.append(np.linalg.solve(chol, x[idx] / sigma[idx, None]))
        whitened_y.append(np.linalg.solve(chol, y[idx] / sigma[idx]))
    a = np.vstack(whitened_x)
    b = np.concatenate(whitened_y)
    coefficients, _, rank, singular = np.linalg.lstsq(a, b, rcond=None)
    if rank != a.shape[1] or singular[-1] <= singular[0] * 1e-10:
        raise ValueError("spectral design is rank deficient or poorly conditioned")
    residual = b - a @ coefficients
    cov = np.linalg.inv(a.T @ a)
    dof = len(points) - len(columns)
    return {
        "n_bins": len(points), "n_products": len(instruments), "dof": dof,
        "rho_assumed": float(rho), "error_floor_ppm": float(error_floor_ppm),
        "chi2": float(residual @ residual),
        "reduced_chi2": float((residual @ residual) / dof),
        "condition_number_design": float(singular[0] / singular[-1]),
        "coefficients_ppm": dict(zip(columns, map(float, coefficients))),
        "conditional_1sigma_ppm": dict(zip(columns, map(float, np.sqrt(np.diag(cov))))),
    }


def leave_one_product_out(points: list[dict], **kwargs) -> dict:
    """Jackknife-style product deletion; NOT independent posterior samples."""
    baseline = fit_spectrum(points, **kwargs)
    outputs = {}
    for name in sorted({p["instrument"] for p in points}):
        retained = [p for p in points if p["instrument"] != name]
        try:
            fitted = fit_spectrum(retained, **kwargs)
            outputs[name] = {
                "n_bins": fitted["n_bins"], "reduced_chi2": fitted["reduced_chi2"],
                "delta_h2o_ppm": fitted["coefficients_ppm"]["h2o"] - baseline["coefficients_ppm"]["h2o"],
                "delta_co2_ppm": fitted["coefficients_ppm"]["co2"] - baseline["coefficients_ppm"]["co2"],
            }
        except ValueError as exc:
            outputs[name] = {"unidentifiable": str(exc)}
    return {"baseline": baseline, "deleted_product": outputs}


def covariance_sweep(points: list[dict], rhos=(0.0, 0.3, 0.6), **kwargs) -> list[dict]:
    """Conditional feature sensitivity to *assumed* within-product covariance."""
    return [fit_spectrum(points, rho=float(rho), **kwargs) for rho in rhos]


def load_benchmark(path: str | Path) -> list[dict]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return payload["points"]
