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

def integrate(x, f, l, n=0, *, sinx=None, cosx=None, kind='auto'):
    """
    Integrate int_0^infty j_l(x) f(x) x^n dx
    """

    assert(x.dtype == 'float64')
    assert(f.dtype == 'float64')

    if not isinstance(n, int):
        raise('n must be an integer: {}' % n)

    if kind == 'trapezoidal':
        _, integ = c._spherical_bessel_integrate('trapezoidal', 1.0, x, f,
                                                 None, None, 0, n, l, 0.0)
        return integ
    elif kind != 'auto':
        raise ValueError('Unknown kind of integral %s' % kind)
    
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
    elif l == 1:
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
    elif l == 2:
        # j2(x) = (3.0/x^3) sin x - 3.0/x^2 cos x - 1/x sin x
        i, integ = c._spherical_bessel_integrate('trapezoidal', 1.0, x, f,
                                                 sinx, cosx, 0, n,
                                                 l, 1.0)
        
        ff = (3.0*x**(n - 3) - 1.0*x**(n - 1))*f
        _, res = c._spherical_bessel_integrate('sin_integ', 1.0, x, ff,
                                               sinx, cosx, i, 0,
                                               0, 0.0)
        integ += res

        ff = -3.0*f*x**(n - 2)
        _, res = c._spherical_bessel_integrate('cos_integ', 1.0, x, ff,
                                               sinx, cosx, i, 0,
                                               0, 0.0)
        integ += res
    elif l == 3:
        # j3(x) = [15/x^4 + 6/x^2 ] sin x
        #         + [-15/x^3 + 1/x] cos x
        i, integ = c._spherical_bessel_integrate('trapezoidal', 1.0, x, f,
                                                 sinx, cosx, 0, n,
                                                 l, 10.0)

        x = x[i:]
        f = f[i:]
        
        ff = (15.0*x**(n - 4) - 6.0*x**(n - 2))*f
        _, res = c._spherical_bessel_integrate('sin_integ', 1.0, x, ff,
                                               sinx[i:], cosx[i:], 0, 0,
                                               0, 0.0)
        integ += res

        ff = (-15.0*x**(n - 3) + x**(n - 1))*f

        _, res = c._spherical_bessel_integrate('cos_integ', 1.0, x, ff,
                                               sinx[i:], cosx[i:], 0, 0,
                                               0, 0.0)
        integ += res
    elif l == 4:
        # j4(x) = [10/x^2 - 105/x^4 ] cos x
        #         + [105/x^5 -45/x^3 + 1/x] sin x
        i, integ = c._spherical_bessel_integrate('trapezoidal', 1.0, x, f,
                                                 sinx, cosx, 0, n,
                                                 l, 10.0)

        x = x[i:]
        f = f[i:]
        
        ff = (105.0*x**(n - 5) - 45.0*x**(n - 3) + 1.0*x**(n - 1))*f
        _, res = c._spherical_bessel_integrate('sin_integ', 1.0, x, ff,
                                               sinx[i:], cosx[i:], 0, 0,
                                               0, 0.0)
        integ += res

        ff = (-105.0*x**(n - 4) + 10.0*x**(n - 2))*f

        _, res = c._spherical_bessel_integrate('cos_integ', 1.0, x, ff,
                                               sinx[i:], cosx[i:], 0, 0,
                                               0, 0.0)
        integ += res
    else:
        raise ValueError('l == %d is not implemented' % l)

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
