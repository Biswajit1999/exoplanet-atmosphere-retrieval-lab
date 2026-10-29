"""Synthetic-only validation; no real atmospheric detection is claimed."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spectral_robustness import (  # noqa: E402
    covariance_sweep, fit_spectrum, leave_one_product_out, templates,
)


def synthetic(seed=4):
    rng = np.random.default_rng(seed)
    points = []
    for name, offset in (("NIRISS", 20500.), ("PRISM", 20800.), ("G395H", 20700.)):
        wavelengths = np.linspace(.7, 5.2, 75)
        h, c, s = templates(wavelengths)
        depth = offset + 15 * np.log(wavelengths / 2.2) + 410 * h + 690 * c + 40 * s
        for w, y in zip(wavelengths, depth + rng.normal(0, 12, len(depth))):
            points.append({"wavelength_um": float(w), "transit_depth_ppm": float(y),
                           "uncertainty_ppm": 12., "instrument": name})
    return points


class SpectralRobustnessTests(unittest.TestCase):
    def test_recovers_injected_contrasts_and_offsets(self):
        result = fit_spectrum(synthetic())
        self.assertAlmostEqual(result["coefficients_ppm"]["h2o"], 410, delta=18)
        self.assertAlmostEqual(result["coefficients_ppm"]["co2"], 690, delta=18)
        self.assertAlmostEqual(result["coefficients_ppm"]["offset:PRISM"] - result["coefficients_ppm"]["offset:NIRISS"], 300, delta=8)
        self.assertTrue(.5 < result["reduced_chi2"] < 1.5)

    def test_noiseless_recovery(self):
        points = synthetic()
        # Restore the exact generative signal rather than fitting a noisy fixture.
        for p in points:
            w = np.array([p["wavelength_um"]]); h, c, s = templates(w)
            offset = {"NIRISS": 20500., "PRISM": 20800., "G395H": 20700.}[p["instrument"]]
            p["transit_depth_ppm"] = float(offset + 15 * np.log(w[0] / 2.2) + 410 * h[0] + 690 * c[0] + 40 * s[0])
        fit = fit_spectrum(points)
        self.assertAlmostEqual(fit["coefficients_ppm"]["h2o"], 410, places=6)
        self.assertLess(fit["chi2"], 1e-15)

    def test_correlation_changes_conditional_uncertainty(self):
        results = covariance_sweep(synthetic(), rhos=(0., .4, .7))
        errors = [r["conditional_1sigma_ppm"]["h2o"] for r in results]
        self.assertGreater(errors[-1], errors[0])
        self.assertTrue(all(np.isfinite(errors)))

    def test_leave_one_product_out_noiseless(self):
        results = leave_one_product_out(synthetic())
        self.assertEqual(len(results["deleted_product"]), 3)
        self.assertTrue(all(abs(v["delta_h2o_ppm"]) < 30 for v in results["deleted_product"].values()))

    def test_rejects_nonphysical_inputs(self):
        p = synthetic()
        with self.assertRaises(ValueError):
            fit_spectrum(p, rho=1.)
        with self.assertRaises(ValueError):
            fit_spectrum(p, error_floor_ppm=-1)
        p[0]["uncertainty_ppm"] = 0
        with self.assertRaises(ValueError):
            fit_spectrum(p)

    def test_rank_deficiency_fails_closed(self):
        p = synthetic()
        for q in p:
            q["wavelength_um"] = 2.0
        with self.assertRaises(ValueError):
            fit_spectrum(p)


if __name__ == "__main__":
    unittest.main()
