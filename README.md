# Exoplanet Atmosphere Retrieval Lab

Transmission spectra, molecular bands and degeneracy diagnostics.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `scripts/validate.js`: no-dependency repository validation.

## Run

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and absence of unfinished scaffold tokens.

## Reference Data

Representative water-band transmission-spectrum anchors for validating molecular-band rendering.

## References

- Seager, S. and Sasselov, D.D., 2000. Theoretical transmission spectra during extrasolar giant planet transits. The Astrophysical Journal, 537(2), pp.916-921.
- Deming, D. et al., 2013. Infrared transmission spectroscopy of the exoplanets HD 209458b and XO-1b using the Wide Field Camera-3 on the Hubble Space Telescope. The Astrophysical Journal, 774(2), p.95.

## Research Quality Upgrade

See [RESEARCH_QUALITY.md](RESEARCH_QUALITY.md) for the validation layer, reference anchors, equations and research boundaries added to this repository.
