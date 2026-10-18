# WASP-39 b conditional template robustness (Python)

This module complements, but does not replace, the JavaScript 51x51 nonnegative-amplitude grid in `inference.js`. It independently solves the **unconstrained** weighted linear regression with the *same Gaussian feature templates*, six per-product offsets, optional log-wavelength slope, and signed SO2 nuisance term. Its coefficients are in ppm of template contrast, **not** molecular abundances, mixing ratios, or detection significances.

The three scientific checks are:

1. **Rank-aware weighted least squares and formal covariance.** Reports template coefficients, nominal 1-sigma errors, design conditioning, chi-square and degrees of freedom. It assumes independent bins; this is not a reconstruction of unavailable JWST inter-bin covariance.
2. **Error-floor sensitivity.** Refit with additional unverified white-noise floors of 0, 50, 100, 200 and 300 ppm. These are stress tests, not calibrated uncertainty inflation.
3. **Leave-one-instrument-out sensitivity.** Refit after excluding each product to identify which spectral regions dominate template coefficients. A missing spectral band can make coefficients unstable or nonidentifiable.

```bash
python -m pip install -r requirements-research.txt
python -m pytest -q tests/test_spectral_robustness.py
python scripts/spectral_robustness.py data/wasp39b_benchmark.json --output results/spectral_robustness.json
```

An independent numerical check on the 355 committed benchmark bins reproduced a baseline unconstrained H2O contrast near 670.8 ppm, CO2 near 979.2 ppm, and reduced chi-square about 8.51. Adding 100 ppm white-noise floor lowered the conditional reduced chi-square to about 3.03; adding 300 ppm lowered it to about 0.59. Removing NIRISS SOSS order 1 shifted the fitted H2O contrast by approximately +997 ppm, demonstrating model/instrument sensitivity. These are not new JWST observations or credible abundance constraints.

Provenance and scientific boundaries follow the pinned Carter/May et al. (2024) uniformly reduced WASP-39 b benchmark spectrum in the existing repository README and `data/wasp39b_benchmark.json`.
