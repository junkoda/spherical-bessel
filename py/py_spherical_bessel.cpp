//#include <iostream> // debug
#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include "sincos_integ.h"
#include "py_spherical_bessel.h"

using namespace std;

//
//
//
static int integrate_trapezoidal(double k,
				 double* r, const size_t r_stride,
				 double* f, const size_t f_stride,
				 const size_t iend,
				 const int l, const int n,
				 const double kr_max,
				 double* const integ);

static double integrate_sin(double k,
			    double* r, const size_t r_stride,
			    double* f, const size_t f_stride,
			    const double sinkr[], const double coskr[],
			    const size_t ibegin, const size_t iend,
			    const int n);

static double integrate_cos(double k,
			    double* r, const size_t r_stride,
			    double* f, const size_t f_stride,
			    const double sinkr[], const double coskr[],
			    const size_t ibegin, const size_t iend,
			    const int n);

class TypeError{};

static void decode_array(const char name[],
			 PyObject* py_obj, Py_buffer* buf,
			 Py_ssize_t len=0, int flags=0)
{
  char msg[128];
  
  if(PyObject_GetBuffer(py_obj, buf, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    throw TypeError();
  
  if(buf->ndim != 1) {
    sprintf(msg, "Expected a 1-dimensional array for %s", name);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(strcmp(buf->format, "d") != 0) {
    sprintf(msg, "Expected an array of double for %s: %s", name, buf->format);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }

  if(len > 0 && buf->shape[0] != len) {
    sprintf(msg, "Expected the length arrays of %d for %s: %d",
	    (int) len, name, (int) buf->shape[0]);
    PyErr_SetString(PyExc_TypeError, msg);
    throw TypeError();
  }
}

//
// Python interface
//
PyObject* py_spherical_bessel_integrate(PyObject* self, PyObject* args)
{
  // _spherical_bessel_integrate(cmd, k, r, f, sin, cos, ibegin, n, l)
  // Integrate: \int_0^\infty d(kr) j_l(kr) (kr)^n f(x)
  //
  // Args:
  //   cmd (string): trapezoidal, sin_integ, cos_integ
  //   k (float)
  //   r (array): r[i] = r_i
  //   f (array): f[i] = f(r_i)
  //   sinx (array): sinx[i] = sin(k r_i), None for trapezoidal
  //   cosx (array): cosx[i] = cos(k r_i), None for trapezoidal
  //   ibegin (int): array index range begin
  //   n (int): n >= 0
  //   l (int): l >= 0 (0 for sin_integ, cos_integ)

  double k;
  int n, l;
int ibegin;
  PyObject *py_r, *py_f, *py_sin, *py_cos;
  char* cmd;
  double kr_max;

  if(!PyArg_ParseTuple(args, "sdOOOOiiid",
		       &cmd, &k,
		       &py_r, &py_f, &py_sin, &py_cos,
		       &ibegin, &n, &l, &kr_max))
    return NULL;
     

  //
  // Decode array information
  //
  Py_buffer r, f, sinx, cosx;
  sinx.buf= cosx.buf= 0;

  try {
    decode_array("r", py_r, &r);
    
    decode_array("f", py_f, &f, r.shape[0]);
  
    if(py_sin != Py_None) 
      decode_array("sinx", py_sin, &sinx, r.shape[0], PyBUF_ANY_CONTIGUOUS);

    if(py_cos != Py_None)
      decode_array("cosx", py_cos, &cosx, r.shape[0], PyBUF_ANY_CONTIGUOUS);
  }
  catch(TypeError) {
    //PyErr_SetString(PyExc_ValueError, msg);
    return NULL;
  }
  

  int i= r.shape[0];
  double integ= 0.0;
  
  if(strcmp(cmd, "trapezoidal") == 0) {
    i= integrate_trapezoidal(k,
			  (double*) r.buf, r.strides[0],
			  (double*) f.buf, f.strides[0],
			  r.shape[0],
			  l, n, kr_max, &integ);
  }
  else if(strcmp(cmd, "sin_integ") == 0) {
    assert(sinx.buf); assert(cosx.buf); assert(l == 0);

    integ = integrate_sin(k,
			  (double*) r.buf, r.strides[0],
			  (double*) f.buf, f.strides[0],
			  (double*) sinx.buf, (double*) cosx.buf,
			  ibegin, r.shape[0], n);
  }
  else if(strcmp(cmd, "cos_integ") == 0) {
    assert(sinx.buf); assert(cosx.buf); assert(l == 0);

    integ = integrate_cos(k,
			  (double*) r.buf, r.strides[0],
			  (double*) f.buf, f.strides[0],
			  (double*) sinx.buf, (double*) cosx.buf,
			  ibegin, r.shape[0], n);
  }
  else {
    char msg[64];
    sprintf(msg, "Unknown command: %.60s", cmd);
    PyErr_SetString(PyExc_ValueError, msg);
    return NULL;
  }
  
  PyBuffer_Release(&r);
  PyBuffer_Release(&f);

  if(py_sin != Py_None)
    PyBuffer_Release(&sinx);
  if(py_cos != Py_None)
    PyBuffer_Release(&cosx);

  
  return Py_BuildValue("id", i, integ);
}

//
// Integration
//

// Condition switching from trapezoidal to piece-wise liear * trigonometric
static inline bool condition(const double kr1, const double kr2,
			     const double kr_max)
{
  //if(kr2 - kr1 > 0.01)
  //if(kr2 > 10.0)
  
  if(kr2 > kr_max)
    return true;
  
  return false;
}


int integrate_trapezoidal(double k,
			  double* r, const size_t r_stride,
			  double* f, const size_t f_stride,
			  const size_t iend,
			  const int l, const int n,
			  const double kr_max,
			  double* const integ)
{
  // NOTE: kr_max is negelected now
  assert(l >= 0);
  // int_0^kr_max (kr)^n j_\ell(kr) f(r)
  double kr1= k*(*r);
  double y1;

  if(l == 0) { // j0
    y1= gsl_sf_bessel_j0(kr1)*pow(kr1, n)*(*f);
      
    for(size_t i=1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= gsl_sf_bessel_j0(kr2)*pow(kr2, n)*(*f);

      *integ += 0.5*(y1 + y2)*(kr2 - kr1);

      if(condition(kr1, kr2, kr_max)) return i;
      
      kr1= kr2;
      y1= y2;
    }
  }
  else if(l == 1) { // j1
    y1= gsl_sf_bessel_j1(kr1)*pow(kr1, n)*(*f);
      
    for(size_t i=1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= gsl_sf_bessel_j1(kr2)*pow(kr2, n)*(*f);

      *integ += 0.5*(y1 + y2)*(kr2 - kr1);

      if(condition(kr1, kr2, kr_max)) return i;

      kr1= kr2;
      y1= y2;      
    }
  }
  else if(l == 2) { // j2
    y1= gsl_sf_bessel_j2(kr1)*pow(kr1, n)*(*f);
      
    for(size_t i=1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= gsl_sf_bessel_j2(kr2)*pow(kr2, n)*(*f);

      *integ += 0.5*(y1 + y2)*(kr2 - kr1);

      if(condition(kr1, kr2, kr_max)) return i;
      
      kr1= kr2;
      y1= y2;      
    }
  }
  else {
    y1= gsl_sf_bessel_jl(l, kr1)*pow(kr1, n)*(*f);
      
    for(size_t i=1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= gsl_sf_bessel_jl(l, kr2)*pow(kr2, n)*(*f);

      *integ += 0.5*(y1 + y2)*(kr2 - kr1);

      if(condition(kr1, kr2, kr_max)) return i;
      
      kr1= kr2;
      y1= y2;      
    }
  }

  return iend;
}

double integrate_sin(double k,
		     double* r, const size_t r_stride,
		     double* f, const size_t f_stride,
		     const double sinkr[], const double coskr[],
		     const size_t ibegin, const size_t iend,
		     const int n)
{
  assert(0 <= n && n < 3);
  //
  // \int (kr)^n dkr sin(kr) f(r)
  //
  //
  r = (double*) ((char*) r + ibegin*r_stride);
  f = (double*) ((char*) f + ibegin*f_stride);

  double kr1= k*(*r);
  double y1= *f;

  double integ= 0.0;
  
  if(n == 0) {    
    for(size_t i=ibegin + 1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= *f;

      // interpolate y = a0 + a1*(kr)

      double a0= (y1*kr2 - y2*kr1)/(kr2 - kr1);
      double a1= (y2 - y1)/(kr2 - kr1);

      integ += a1*(  sin_integ1(kr2, sinkr[i],   coskr[i])
		   - sin_integ1(kr1, sinkr[i-1], coskr[i-1]))
	     + a0*(  sin_integ0(kr2, sinkr[i],   coskr[i])
		   - sin_integ0(kr1, sinkr[i-1], coskr[i-1]));
      
      kr1= kr2;
      y1= y2;
    }
  }
  else if(n == 1) {
    for(size_t i=ibegin + 1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= *f;

      // interpolate y = a0 + a1*(kr)
      double a0= (y1*kr2 - y2*kr1)/(kr2 - kr1);
      double a1= (y2 - y1)/(kr2 - kr1);

      // y = 
      integ += a1*(  sin_integ2(kr2, sinkr[i],   coskr[i])
		   - sin_integ2(kr1, sinkr[i-1], coskr[i-1]))
	     + a0*(  sin_integ1(kr2, sinkr[i],   coskr[i])
		   - sin_integ1(kr1, sinkr[i-1], coskr[i-1]));

      kr1= kr2;
      y1= y2;
    }
  }
  else if(n == 2) {
    for(size_t i=ibegin + 1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= *f;

      // interpolate y = a0 + a1*(kr)
      double a0= (y1*kr2 - y2*kr1)/(kr2 - kr1);
      double a1= (y2 - y1)/(kr2 - kr1);

      // y = 
      integ += a1*(  sin_integ3(kr2, sinkr[i],   coskr[i])
		   - sin_integ3(kr1, sinkr[i-1], coskr[i-1]))
	     + a0*(  sin_integ2(kr2, sinkr[i],   coskr[i])
		   - sin_integ2(kr1, sinkr[i-1], coskr[i-1]));

      kr1= kr2;
      y1= y2;
    }
  }

  return integ;
}


double integrate_cos(double k,
		     double* r, const size_t r_stride,
		     double* f, const size_t f_stride,
		     const double sinkr[], const double coskr[],
		     const size_t ibegin, const size_t iend,
		     const int n)
{
  assert(0 <= n && n < 3);
  //
  // \int r^n dr j_l(kr) f(r)
  //
  //
  r = (double*) ((char*) r + ibegin*r_stride);
  f = (double*) ((char*) f + ibegin*f_stride);

  double kr1= k*(*r);
  double y1= *f;

  double integ= 0.0;
  
  if(n == 0) {
    for(size_t i=ibegin + 1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= *f;

      // interpolate y = a0 + a1*(kr)
      double a0= (y1*kr2 - y2*kr1)/(kr2 - kr1);
      double a1= (y2 - y1)/(kr2 - kr1);

      // y = 
      integ += a1*(  cos_integ1(kr2, sinkr[i],   coskr[i])
		   - cos_integ1(kr1, sinkr[i-1], coskr[i-1]))
	     + a0*(  cos_integ0(kr2, sinkr[i],   coskr[i])
		   - cos_integ0(kr1, sinkr[i-1], coskr[i-1]));

      kr1= kr2;
      y1= y2;
    }
  }
  else if(n == 1) {
    for(size_t i=ibegin + 1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= *f;

      // interpolate y = a0 + a1*(kr)
      double a0= (y1*kr2 - y2*kr1)/(kr2 - kr1);
      double a1= (y2 - y1)/(kr2 - kr1);

      // y = 
      integ += a1*(  cos_integ2(kr2, sinkr[i],   coskr[i])
		   - cos_integ2(kr1, sinkr[i-1], coskr[i-1]))
	     + a0*(  cos_integ1(kr2, sinkr[i],   coskr[i])
		   - cos_integ1(kr1, sinkr[i-1], coskr[i-1]));

      kr1= kr2;
      y1= y2;
    }
  }
  else if(n == 2) {
    for(size_t i=ibegin + 1; i<iend; ++i) {
      r = (double*) ((char*) r + r_stride);
      f = (double*) ((char*) f + f_stride);
      double kr2= k*(*r);
      double y2= *f;

      // interpolate y = a0 + a1*(kr)
      double a0= (y1*kr2 - y2*kr1)/(kr2 - kr1);
      double a1= (y2 - y1)/(kr2 - kr1);

      // y = 
      integ += a1*(  cos_integ3(kr2, sinkr[i],   coskr[i])
		   - cos_integ3(kr1, sinkr[i-1], coskr[i-1]))
	     + a0*(  cos_integ2(kr2, sinkr[i],   coskr[i])
		   - cos_integ2(kr1, sinkr[i-1], coskr[i-1]));

      kr1= kr2;
      y1= y2;
    }
  }

  return integ;
}


