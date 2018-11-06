#include <cmath>
#include "py_spherical_bessel.h"

//
//
//
//py_spherical_bessel_module_init();

//
// Python interface
//
PyObject* py_spherical_bessel_integrate(PyObject* self, PyObject* args)
{
  // _spherical_bessel_integrate(k, r, f, n, l)
  // Integrate: \int_0^\infty d(kr) j_l(kr) (kr)^n f(x)
  //
  // Args:
  //   k (float)
  //   r (array): r[i] = r_i
  //   f (array): f[i] = f(r_i)
  //   n (int): n >= 0
  //   l (int): l >= 0

  double k;
  int n, l;
  PyObject *py_r, *py_f;

  if(!PyArg_ParseTuple(args, "dOOii",
		       &k, &py_r, &py_f, &n, &l))
    return NULL;

    //
  // Decode array information
  //

  // r array
  Py_buffer r;
  if(PyObject_GetBuffer(py_r, &r, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;
  
  if(r.ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected a 1-dimensional array for r");
    return NULL;
  }

  if(!strcmp(r.format, "d")) {
    PyErr_SetString(PyExc_TypeError, "Expected an array of double for r");    
  }

  // f array
  Py_buffer f;
  if(PyObject_GetBuffer(py_f, &f, PyBUF_FORMAT | PyBUF_FULL_RO) == -1)
    return NULL;

  if(r.ndim != 1) {
    PyErr_SetString(PyExc_TypeError, "Expected a 1-dimensional array f");
    return NULL;
  }

  if(!strcmp(r.format, "f")) {
    PyErr_SetString(PyExc_TypeError, "Expected an array of double for f");    
  }


  double integ= 0.0;
  
  return Py_BuildValue("d", integ);
}

  

