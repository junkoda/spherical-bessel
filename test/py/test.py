import numpy as np
import math
import unittest

import spherical_bessel as s

class TestSphericalBessel(unittest.TestCase):
    """
    Test spherical_bessel library
    """

    def test_j0(self):
        """
        Test S j0(x) dx = pi/2
        """

        x = 10**np.linspace(-4, 5, 1001)
        f = np.ones_like(x)
        
        self.assertAlmostEqual(s.integrate(x, f, 0), math.pi/2, places=4)


if __name__ == "__main__":
    unittest.main()
