//
// Information for Python
//
#include <iostream>
#include "Python.h"

#include "py_spherical_bessel.h"

using namespace std;

static PyMethodDef methods[] = {
  {"_spherical_bessel_integrate",
   py_spherical_bessel_integrate, METH_VARARGS,
   "integrate S j_l (kr)^2 j(kr) dr"},
  
  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef module = {
  PyModuleDef_HEAD_INIT,
  "_spherical_bessel", // name of this module
  "Spherical bessel integral module", // Doc String
  -1,
  methods
};

PyMODINIT_FUNC
PyInit__spherical_bessel(void) {
  //PyInit__spherical_bessel(void) {
  //py_spherical_bessel_module_init();
  
  return PyModule_Create(&module);
}
