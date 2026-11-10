# Conditional spectral-inference diagnostics (2026-11-10)

Standalone Python numerical cross-checks of the browser's reduced-order WASP-39 b model. Requires NumPy; run `cd research && python -m unittest -v test_robust_inference.py`.

- `gls_fit`: Cholesky whitening followed by SVD least squares, with coefficient covariance and chi-square; rejects non-SPD covariance and non-identifiable design.
- `ar1_covariance`: a heteroscedastic, **assumed** adjacent-bin AR(1) correlation model. Bin ordering must be physically meaningful.
- `profile_template`: template amplitude and conditional uncertainty after profiling nuisance offsets/slopes.

Tests use only synthetic arrays, not the benchmark WASP-39 b spectrum. Correlation is not measured, the templates are not molecular opacities, and this is not an abundance retrieval or detection significance. Full benchmark integration and covariance provenance remain future work. AI-assisted implementation; Biswajit Jana is the requested Git author.
