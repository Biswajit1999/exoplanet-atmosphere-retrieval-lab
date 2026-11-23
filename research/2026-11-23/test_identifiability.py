import unittest
import numpy as np
from identifiability import singular_diagnostics, linear_estimate, gaussian_coverage

class TestIdentifiability(unittest.TestCase):
    def setUp(self):
        t = np.linspace(-1, 1, 40)
        self.x = np.column_stack([np.ones(len(t)), t])
        self.s = np.full(len(t), 0.1)
        self.truth = np.array([2., -0.5])

    def test_fit_and_covariance(self):
        beta, cov = linear_estimate(self.x, self.x @ self.truth, self.s)
        np.testing.assert_allclose(beta, self.truth, atol=1e-12)
        np.testing.assert_allclose(cov, np.linalg.inv((self.x/self.s[:,None]).T @ (self.x/self.s[:,None])), atol=1e-12)

    def test_conditioning(self):
        singular, condition = singular_diagnostics(self.x, self.s)
        self.assertGreater(singular[-1], 0)
        self.assertLess(condition, 3)

    def test_rank_failure(self):
        with self.assertRaises(ValueError):
            linear_estimate(np.ones((40, 2)), np.ones(40), self.s)

    def test_coverage_calibration(self):
        cover, spread = gaussian_coverage(self.x, self.s, self.truth, repetitions=800, seed=17)
        np.testing.assert_array_less(np.abs(cover - 0.95), np.full(2, 0.045))
        np.testing.assert_array_less(np.abs(spread - 1), np.full(2, 0.12))

    def test_misspecified_noise(self):
        cover, spread = gaussian_coverage(self.x, self.s, self.truth, repetitions=800, seed=17, noise_scale=2)
        self.assertTrue(np.all(cover < 0.75))
        self.assertTrue(np.all(spread > 1.7))

if __name__ == "__main__":
    unittest.main()
