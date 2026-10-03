# WASP-39 b Spectral Inference Laboratory

[![Dataset DOI](https://img.shields.io/badge/data-10.5281%2Fzenodo.10161743-0b665f)](https://doi.org/10.5281/zenodo.10161743)
[![Publication DOI](https://img.shields.io/badge/paper-10.1038%2Fs41550--024--02292--x-a65c24)](https://doi.org/10.1038/s41550-024-02292-x)

A browser-based, reduced-order likelihood analysis of the uniformly reduced
JWST transmission spectrum of WASP-39 b. The laboratory is designed to show
both what a simple molecular-feature model extracts and where that model
fails. It is not a line-by-line radiative-transfer retrieval and does not
convert feature amplitude into molecular abundance.

## Scientific question

How much H2O-like and CO2-like spectral contrast is preferred by the public
0.52–5.35 micrometre benchmark spectrum after allowing every detector/order
product its own transit-depth offset, and how sensitive is that result to a
spectral slope, template width, and unmodelled error floor?

The default analysis answers this conditional question with all six products,
published marginal uncertainties, no extra error floor, and a log-wavelength
slope. It also reports a goodness-of-fit check before any interpretation.

## Data

The source is Carter, May et al. (2024), *A benchmark JWST near-infrared
spectrum for the exoplanet WASP-39 b*, Nature Astronomy 8, 1008–1019. The
complete upstream Zenodo archive is 1.7 GB and has MD5
`f617fa28908aa35d334062f1c06bc24d`. This repository pins the six small source
tables used by the browser:

| Product | Selection | Included bins |
|---|---:|---:|
| NIRISS/SOSS order 1 | R=100 | 108 |
| NIRISS/SOSS order 2 | R=100 | 27 |
| NIRCam/F322W2 | R=100 | 49 |
| NIRSpec/G395H NRS1 | R=100 | 30 |
| NIRSpec/G395H NRS2 | R=100 | 29 |
| NIRSpec/PRISM | native `bins_scale1` | 112 |

Thirty-five PRISM bins with a non-zero upstream `satflag` are excluded. Every
source CSV is checked against the SHA-256 digest embedded in
`data/wasp39b_benchmark.json`.

The source tables report radius ratio. The committed browser dataset applies

```text
depth_ppm = (Rp/Rs)^2 * 10^6
sigma_depth_ppm = 2 * (Rp/Rs) * mean(sigma_low, sigma_high) * 10^6
```

The uncertainty conversion is a first-order propagation. It symmetrizes the
reported asymmetric radius-ratio errors and does not reconstruct inter-bin
covariance.

## Declared model

For spectral bin `i` from product `g`, the laboratory evaluates

```text
y_i = b_g
    + s log(lambda_i / 2.2 micrometre)
    + A_H2O T_H2O(lambda_i)
    + A_CO2 T_CO2(lambda_i)
    + a_SO2 T_SO2(lambda_i)
    + epsilon_i

epsilon_i ~ Normal(0, sigma_i^2 + f^2)
```

`b_g` is a separately profiled offset for each detector/order product. `s` is
an optional continuum slope. `a_SO2` is a signed nuisance/decorrelation
coefficient, not an SO2 abundance. `f` is the user-controlled additional
error floor. The H2O template is a normalized sum of Gaussian bands near
1.15, 1.40, 1.90, and 2.70 micrometres; CO2 and SO2 use simple features near
4.30 and 4.05 micrometres. A width-scale control moves all template widths
together for sensitivity analysis.

At every point on a 51 x 51 `(A_H2O, A_CO2)` grid, the code solves the weighted
normal equations for all nuisance coefficients. The grid uses independent
flat priors on each non-negative amplitude from 0 to 1400 ppm. It reports
normalized likelihood weights and equal-tail 16th, 50th, and 84th percentiles.
The displayed 2-D contours are likelihood-ratio guides, not calibrated
coverage regions under model misspecification.

## Default result and its interpretation

The committed 355-bin configuration gives approximately:

- H2O-template amplitude: 672 ppm, conditional grid interval [644, 700] ppm;
- CO2-template amplitude: 980 ppm, conditional grid interval [952, 1008] ppm;
- residual RMS: 250.4 ppm;
- chi-squared / dof: 2935.8 / 345, or reduced chi-squared 8.51.

The last line dominates the interpretation. Under the published marginal
errors, this template family is not an adequate generative model of the
benchmark spectrum. Consequently, its narrow conditional intervals and its
large BIC contrast against a nuisance model with H2O and CO2 amplitudes fixed
to zero must **not** be reported as abundance constraints or molecular
detection significances. The result is an
auditable demonstration of feature contrast and model inadequacy.

This is also why the interface exposes an error-floor stress test. An added
error floor can make reduced chi-squared approach unity, but doing so is a
sensitivity experiment, not evidence that the extra variance is physically
correct. A research-grade abundance retrieval would require wavelength-
dependent opacity, radiative transfer, chemically meaningful priors,
instrument systematics, and preferably a released covariance model.

## Reproduce

The website has no runtime dependencies. Python is required only to rebuild
the compact JSON from the pinned CSV products.

```bash
python scripts/build_dataset.py
npm run check
python -m http.server 8080
```

Open `http://localhost:8080`. The numerical grid runs in a Web Worker so the
controls and accessibility tree remain responsive.

`npm run check` performs syntax checks, validates all six product hashes and
355 finite bins, checks the public inference boundaries, and runs four
numerical tests. The tests recover injected amplitudes and offsets, verify
error-floor behaviour, and ensure undersized data selections fail closed.

## Repository map

```text
data/source_products/       six published spectrum tables
data/wasp39b_benchmark.json browser-ready data plus provenance
scripts/build_dataset.py    deterministic radius-ratio conversion
scripts/validate.js         integrity and public-boundary checks
inference.js                pure likelihood and profiling implementation
physicsWorker.js            worker boundary
app.js                      accessible charts and interaction
tests/inference.test.js     numerical regression tests
```

## Inference boundaries

- The selected CSV products do not supply a joint inter-bin covariance matrix;
  the likelihood is diagonal with an optional scalar error floor.
- Overlapping instruments are not independent demonstrations of the same
  molecule, and per-product offsets remove absolute-radius information.
- Gaussian feature templates are not opacity cross-sections and their
  amplitudes are not volume-mixing ratios, metallicity, or C/O.
- The grid boundaries and flat priors affect conditional intervals.
- The browser fit is not independent validation of the Carter et al. reduction;
  it reuses their released spectrum.
- Published molecular claims should be taken from full atmospheric analyses,
  not from this reduced-order model.

## References

- Carter, A. L., May, E. M. et al. (2024). A benchmark JWST near-infrared
  spectrum for the exoplanet WASP-39 b. *Nature Astronomy*, 8, 1008–1019.
  <https://doi.org/10.1038/s41550-024-02292-x>
- JWST Transiting Exoplanet Community Early Release Science Team (2023).
  Identification of carbon dioxide in an exoplanet atmosphere. *Nature*, 614,
  649–652. <https://doi.org/10.1038/s41586-022-05269-w>
- Alderson, L. et al. (2023). Early Release Science of the exoplanet WASP-39b
  with JWST NIRSpec G395H. *Nature*, 614, 664–669.
  <https://doi.org/10.1038/s41586-022-05591-3>
- Feinstein, A. D. et al. (2023). Early Release Science of the exoplanet
  WASP-39b with JWST NIRISS. *Nature*, 614, 670–675.
  <https://doi.org/10.1038/s41586-022-05674-1>

## Licence

Code is MIT licensed. The published data retain their upstream attribution
and licensing; cite the dataset and paper when reusing the spectrum.
