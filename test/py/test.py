import numpy as np
import math
import unittest

import spherical_bessel as s

class TestSphericalBessel(unittest.TestCase):
    """
    Test spherical_bessel library
    """

    def test_sin_integ0(self):
        x = np.array([0.5, 1.5])
        f = np.array([1.0, 1.0])

        test = s.sin_integ(x, f, 0)
        self.assertAlmostEqual(test, 0.806845, places=4)

    def test_sin_integ1(self):
        x = np.array([0.5, 1.5])
        f = np.array([1.0, 1.0])

        test = s.sin_integ(x, f, 1)
        self.assertAlmostEqual(test, 0.850755, places=4)

    def test_sin_integ2(self):
        x = np.array([0.5, 1.5])
        f = np.array([1.0, 1.0])

        test = s.sin_integ(x, f, 2)
        self.assertAlmostEqual(test, 0.959606, places=4)

    def test_cos_integ0(self):
        x = np.array([0.5, 1.5])
        f = np.array([1.0, 1.0])

        test = s.cos_integ(x, f, 0)
        self.assertAlmostEqual(test, 0.518069, places=4)

    def test_cos_integ1(self):
        x = np.array([0.5, 1.5])
        f = np.array([1.0, 1.0])

        test = s.cos_integ(x, f, 1)
        self.assertAlmostEqual(test, 0.449684, places=4)

    def test_cos_integ2(self):
        x = np.array([0.5, 1.5])
        f = np.array([1.0, 1.0])

        test = s.cos_integ(x, f, 2)
        self.assertAlmostEqual(test, 0.422997, places=4)

    def test_j0(self):
        """
        Test int j0(x) dx = pi/2
        """

        x = 10**np.linspace(-4, 5, 1001)
        f = np.ones_like(x)
        
        self.assertAlmostEqual(s.integrate(x, f, 0), math.pi/2, places=4)

    def test_j1(self):
        """
        Test int j1(x) dx = 1
        """

        x = 10**np.linspace(-4, 5, 2001)
        f = np.ones_like(x)
        
        self.assertAlmostEqual(s.integrate(x, f, 1), 1.0, places=4)

    def test_j2(self):
        """
        Test int j2(x) dx = 0.785398
        """

        x = 10**np.linspace(-4, 5, 3001)
        f = np.ones_like(x)
        
        self.assertAlmostEqual(s.integrate(x, f, 2), 0.785398, places=4)

    def test_j3(self):
        """
        Test int_1.0e-4^10 j3(x) dx = 0.573031
        """

        x = 10**np.linspace(-4, 5, 4001) 
        f = np.ones_like(x)
        
        self.assertAlmostEqual(s.integrate(x, f, 3), 2.0/3.0, places=4)

    def test_j4(self):
        """
        Test int j2(x) dx = 0.785398
        """

        x = 10**np.linspace(-4, 1.0, 2001)
        f = np.ones_like(x)
        
        self.assertAlmostEqual(s.integrate(x, f, 4), 0.602723, places=4)


if __name__ == "__main__":
    unittest.main()
