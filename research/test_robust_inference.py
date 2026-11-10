import unittest
import numpy as np
from robust_inference import gls_fit, ar1_covariance, profile_template

class TestInference(unittest.TestCase):
    def test_exact_recovery(self):
        x = np.column_stack((np.ones(15), np.linspace(-1, 1, 15)))
        beta, cov, chi, dof = gls_fit(x, x @ np.array([2., -0.5]), np.eye(15))
        np.testing.assert_allclose(beta, [2, -.5], atol=1e-12)
        self.assertLess(chi, 1e-25)
        self.assertEqual(dof, 13)
        self.assertTrue(np.all(np.linalg.eigvalsh(cov)>0))
    def test_ar1(self):
        c = ar1_covariance(np.ones(9)*2, .6)
        self.assertAlmostEqual(c[0, 1], 2.4)
        self.assertTrue(np.all(np.linalg.eigvalsh(c)>0))
    def test_profile(self):
        t = np.linspace(-1,1,30)**2
        a, err, chi, dof = profile_template(np.ones((30,1)), t, 3 + 8*t, np.eye(30))
        self.assertAlmostEqual(a,8,places=11)
        self.assertGreater(err,0)
        self.assertLess(chi,1e-24)
    def test_bad_covariance(self):
        with self.assertRaises(ValueError):
            gls_fit(np.ones((2,1)),np.ones(2),np.ones((2,2)))
    def test_rank(self):
        with self.assertRaises(ValueError):
            gls_fit(np.ones((5,2)),np.ones(5),np.eye(5))
if __name__ == '__main__': unittest.main()
