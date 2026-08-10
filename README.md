# Exoplanet Atmosphere Retrieval Lab

Transmission spectra, molecular bands and degeneracy diagnostics.

Created and maintained by Biswajit Jana.

This is a zero-build, browser-only laboratory for exploring how an exoplanet's
transmission spectrum responds to atmospheric parameters (temperature,
metallicity, cloud-top pressure, base transit depth). It overlays a small,
citable bundle of reference transmission-spectrum anchors on top of an
interactive forward model so you can see, by eye, how model choices trade off
against real measurement uncertainty and degeneracy.

**Scope note:** this repository implements a parametric *forward model* (you
choose parameters, it draws the spectrum) plus a toy 2-D response surface. It
does **not** perform a statistical atmospheric retrieval (no MCMC/nested
sampling, no likelihood evaluated against data, no posterior). The heatmap
panel is a fast qualitative proxy for how sensitive the spectrum is to two
parameters at once, not a fitted probability surface. See "How it works"
below for exactly what the code computes.

## Scientific Background

### Transmission spectroscopy during transit

When a planet transits its host star, starlight grazes the planet's
terminator and passes through its atmosphere before reaching the observer.
Atoms and molecules in the atmosphere absorb specific wavelengths, so the
apparent planet radius — and therefore the transit depth — varies with
wavelength:

```
depth(lambda) = (R_p(lambda) / R_star)^2
```

At wavelengths coincident with a strong molecular absorption feature, the
planet appears larger (light is blocked higher in the atmosphere), so the
transit is deeper. This differential signal, typically tens to hundreds of
parts per million for a hot Jupiter around a Sun-like star, is the
observable that transmission spectroscopy extracts (Seager & Sasselov 2000).

### Atmospheric scale height

The size of the wavelength-dependent signal scales with the atmospheric
pressure scale height, the characteristic altitude over which pressure falls
by 1/e:

```
H = k_B * T / (mu * g)
```

where `k_B` is Boltzmann's constant, `T` the atmospheric temperature, `mu`
the mean molecular mass, and `g` the planet's surface gravity. Hotter,
lower-gravity, lower-molecular-weight atmospheres (hot, low-mass, hydrogen-
dominated planets) produce the largest spectral features and are the most
favorable retrieval targets — a large part of why hot Jupiters and sub-
Neptunes dominate the current atmospheric-characterization sample.

### Molecular absorption features

Different molecules imprint absorption bands at characteristic
near/mid-infrared wavelengths:

- **H2O** — bands near 1.15, 1.4, 1.9 and 2.7 micron; the workhorse feature
  for hot-Jupiter transmission spectroscopy.
- **CH4** — bands near 1.65, 2.3 and 3.3 micron, prominent in cooler,
  higher-C/O atmospheres.
- **CO2** — a strong band near 4.3 micron, a key JWST/NIRSpec diagnostic
  (e.g. the WASP-39b CO2 detection).
- **Na / K** — narrow optical resonance lines near 0.589 and 0.77 micron
  that probe the upper atmosphere and cloud/haze coverage.

Clouds and hazes mute all of these features by raising the effective
opaque surface (the "cloud-top pressure"), flattening the spectrum — a
classic source of degeneracy between "clear atmosphere, low abundance" and
"cloudy atmosphere, high abundance" solutions.

### Retrieval: forward modeling vs. inversion

A **forward model** takes atmospheric parameters (temperature, composition,
cloud deck, planet/star radii) and predicts a spectrum. An **atmospheric
retrieval** runs the inverse problem: given an observed spectrum with error
bars, it searches parameter space (typically via MCMC or nested sampling,
e.g. `PyMultiNest`, `MultiNest`, `emcee`) for the posterior distribution of
atmospheric parameters that best reproduce the data, explicitly quantifying
degeneracies such as the temperature-metallicity-cloud trade-off
(Madhusudhan 2018 review). This repository implements only the forward-model
half of that pipeline, interactively, in the browser — it is a teaching and
intuition-building tool for how the forward model behaves, not a retrieval
engine.

### Relevance to JWST

JWST's NIRSpec (0.6-5.3 micron) and NIRISS SOSS (0.6-2.8 micron) modes have
delivered the highest-precision transmission spectra to date, resolving H2O,
CO2, SO2 and other features in planets such as WASP-39b and WASP-96b at
precisions well below 100 ppm — a regime where the reference anchors and
sliders in this lab are directly comparable in scale to real data.

## How It Works

The app is a single-page, dependency-free static site:

1. `index.html` renders the mission-control layout; `styles.css` styles it.
2. `app.js` builds the parameter sliders from the `LAB.controls` definition
   (temperature, metallicity scale, cloud-top pressure proxy, base radius
   depth), loads `data/reference.json` (the fixed reference anchors, drawn
   as amber points), and on every slider change posts the current parameters
   to `physicsWorker.js`.
3. `physicsWorker.js` runs off the UI thread. Its `atmosphere(p)` function
   computes:
   - a scale-height proxy `scale = 90 * sqrt(T / 1200) * log10(1 + metallicity)`,
   - a cloud transmission factor `cloud = 1 - cloudTop`,
   - a transit-depth spectrum built from four Gaussian absorption bands
     (approximating H2O at 0.95/1.4 micron, CO at 2.3 micron, CO2 at
     4.3 micron) added on top of the base radius depth, each scaled by
     `cloud * scale`,
   - a small metrics block (`water_band_ppm`, `scale_height_proxy`,
     `cloud_transmission`),
   - and a 72x72 normalized response heatmap (labeled "retrieval
     likelihood" in the UI) that is a smooth analytic function of the two
     heatmap axes, used only as a qualitative sensitivity visualization —
     it is not derived from a fit to the reference data.
4. `app.js` receives the worker's result, draws the model spectrum and the
   reference anchors on a shared canvas (`drawSeries`), renders the heatmap
   (`drawHeatmap`), and updates the telemetry panel.
5. `research-overlay.js` adds a small non-blocking status panel confirming
   the reference-anchor and equation count; it has no effect on the model.
6. `scripts/validate.js` (via `npm run check`) is a dependency-free CI-style
   check: it verifies required files exist, the reference JSON has at least
   four finite anchor points, the README cites every reference required by
   `data/reference.json`, and that no unfinished scaffold tokens remain.
   `scripts/validate_repository.mjs` runs a similar check against
   `data/research-reference.json` for the separate research-quality layer
   described in `RESEARCH_QUALITY.md`.

## Usage

```bash
python -m http.server 8080
```

Open `http://localhost:8080` and drag the sliders. Each change triggers a
new worker run (see the `run NNN` counter); the cyan curve is the model
spectrum, the amber points are the fixed reference anchors, and the right
panel shows the qualitative two-parameter response field.

## Validate

```bash
npm run check              # repository/reference/citation checks
npm run validate:research  # research-quality layer checks
```

## Math Appendix

Transit depth as a function of wavelength (idealized):

```
depth(lambda) = (R_p(lambda) / R_star)^2
```

Pressure scale height:

```
H = k_B * T / (mu * g)
```

Model spectrum actually computed by `physicsWorker.js` (`atmosphere(p)`):

```
scale  = 90 * sqrt(T / 1200) * log10(1 + Z)
cloud  = 1 - cloudTop
band(c, w, a, lambda) = a * exp(-0.5 * ((lambda - c) / w)^2)

depth(lambda) = R0
              + cloud * scale * [        band(1.40, 0.12, 1, lambda)
                                  + 0.75 * band(2.30, 0.18, 1, lambda)
                                  + 1.10 * band(4.30, 0.25, 1, lambda)
                                  + 0.35 * band(0.95, 0.08, 1, lambda) ]
```

where `T` is temperature [K], `Z` is the metallicity slider (dimensionless
scale factor, not [Fe/H] dex), `cloudTop` in [0, 1] is the cloud-top
pressure proxy, and `R0` is the base radius depth in ppm. This is a
deliberately simplified, illustrative band model — it captures the correct
qualitative trends (hotter/higher-metallicity atmospheres show taller
bands, higher clouds mute them) but is not a line-by-line or correlated-k
radiative-transfer calculation.

## Reference Data

`data/reference.json` holds representative water-band transmission-spectrum
anchors (HST/JWST-style hot-Jupiter values) used to visually validate that
the model's molecular-band rendering falls in a physically reasonable range.
`data/research-reference.json` holds a separate, smaller benchmark set used
only by the repository's research-quality validation script.

## References

- Seager, S. and Sasselov, D.D., 2000. Theoretical transmission spectra during extrasolar giant planet transits. The Astrophysical Journal, 537(2), pp.916-921.
- Deming, D. et al., 2013. Infrared transmission spectroscopy of the exoplanets HD 209458b and XO-1b using the Wide Field Camera-3 on the Hubble Space Telescope. The Astrophysical Journal, 774(2), p.95.
- Madhusudhan, N., 2019. Exoplanetary atmospheres: key insights, challenges, and prospects. Annual Review of Astronomy and Astrophysics, 57, pp.617-663.
- JWST Transiting Exoplanet Community Early Release Science Team, 2023. Identification of carbon dioxide in an exoplanet atmosphere. Nature, 614(7949), pp.649-652.
- Wilson, G. et al., 2017. Good enough practices in scientific computing. PLOS Computational Biology, 13(6), p.e1005510.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors, equations and research boundaries added to this repository.
