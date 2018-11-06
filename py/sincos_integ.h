#ifndef SINCOS_INTEG_H
#define SINCOS_INTEG_H 1

static inline double sin_integ0(const double x, const double sinx, const double cosx)
{
  // \int sin(x) dx = sin_integ0(x)
  return -cosx;
}

static inline double sin_integ1(const double x, const double sinx, const double cosx)
{
  // \int x sin(x) dx = sin_integ1(x)
  return sinx - x*cosx;
}

static inline double sin_integ2(const double x, const double sinx, const double cosx)
{
  // \int x^2 sin(x) dx = sin_integ2(x)
  return (2.0 - x*x)*cosx + 2.0*x*sinx;
}

static inline double sin_integ3(const double x, const double sinx, const double cosx)
{
  // \int x^3 sin(x) dx = sin_integ3(x)
  return 3.0*(x*x - 2.0)*sinx - x*(x*x - 6.0)*cosx;
}

static inline double cos_integ0(const double x, const double sinx, const double cosx)
{
  // \int cos(x) dx
  return sinx;
}

static inline double cos_integ1(const double x, const double sinx, const double cosx)
{
  // \int x cos(x) dx
  return x*sinx + cosx;
}

static inline double cos_integ2(const double x, const double sinx, const double cosx)
{
  // \int x^2 cos(x) dx
  return (x*x - 2.0)*sinx + 2.0*x*cosx;
}

static inline double cos_integ3(const double x, const double sinx, const double cosx)
{
  // \int x^3 cos(x) dx
  return (3.0*x*x - 6.0)*cosx + (x*x*x - 6.0*x)*sinx;
}





#endif
