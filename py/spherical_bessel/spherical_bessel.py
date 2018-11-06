import spherical_bessel._spherical_bessel as c
import numpy as np

def sin_integ(x, f, n=0, *, sinx=None, cosx=None):
    if sinx is None: sinx = np.sin(x)
    if cosx is None: cosx = np.cos(x)

    _, integ = c._spherical_bessel_integrate('sin_integ', 1.0, x, f,
                                             sinx, cosx, 0, n,
                                             0, 0.0)
    return integ

def cos_integ(x, f, n=0, *, sinx=None, cosx=None):
    if sinx is None: sinx = np.sin(x)
    if cosx is None: cosx = np.cos(x)

    _, integ = c._spherical_bessel_integrate('cos_integ', 1.0, x, f,
                                             sinx, cosx, 0, n,
                                             0, 0.0)
    return integ

def integrate(x, f, l, n=0, *, sinx=None, cosx=None):
    """
    Integrate int_0^infty j_l(x) f(x) x^n dx
    """

    assert(x.dtype == 'float64')
    assert(f.dtype == 'float64')

    if not isinstance(n, int):
        raise('n must be an integer: {}' % n)
    
    i = 0
    integ = 0.0
    
    if sinx is None:
        sinx = np.sin(x)
    if cosx is None:
        cosx = np.cos(x)

    if l == 0:
        # j0(x) = sin(x)/x
        if 1 <= n <= 2:
            i, integ = c._spherical_bessel_integrate('sin_integ', 1.0, x, f,
                                                     sinx, cosx, 0, n - 1,
                                                     0, 0.0)

        else:
            ff = f*x**(n - 1)
            i, integ = c._spherical_bessel_integrate('sin_integ', 1.0, x, ff,
                                                     sinx, cosx, 0, 0,
                                                     0, 0.0)            

    if l == 1:
        print('j1 integral, n=%d' % n)
        # j1(x) = (sin x)/(x^2) - (cos x)/x
        if 2 <= n <= 3:
            _, integ = c._spherical_bessel_integrate('sin_integ', 1.0, x, f,
                                                     sinx, cosx, 0, n - 2,
                                                     0, 0.0)
            _, res = c._spherical_bessel_integrate('sin_integ', 1.0, x, f,
                                                   sinx, cosx, 0, n - 1,
                                                   0, 0.0)
            integ -= res
        elif n < 2:
            i, integ = c._spherical_bessel_integrate('trapezoidal', 1.0, x, f,
                                                     sinx, cosx, 0, n,
                                                     l, 1.0)

            ff = f*x**(n - 2)
            _, res = c._spherical_bessel_integrate('sin_integ', 1.0, x, ff,
                                                   sinx, cosx, i, 0,
                                                   0, 0.0)
            integ += res

            ff = f*x**(n - 1)
            _, res = c._spherical_bessel_integrate('cos_integ', 1.0, x, ff,
                                                   sinx, cosx, i, 0,
                                                   0, 0.0)
            integ -= res
            

            
    return integ

def integrate_dx(k, r, f, l, n=0):
    """
    Integrate 4pi (-i)^l int_0^infity r^n dr j_l(kr) f(r)

    Args:
      k (float): real number k
      r (array): r[i] = r_i
      f (array): f[i] = f(r_i)
      l (int):   l in spherical Bessel function j_l(kr)
      n (int):   factor in r^n
    """
    
    return 0.0
