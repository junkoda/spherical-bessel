#
# $ make
#

from distutils.core import setup, Extension
#import numpy as np
import os


# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++.
# https://stackoverflow.com/questions/8106258
import distutils.sysconfig
cfg_vars = distutils.sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")

# directories for include -I(idir)
#idirs = os.environ["IDIRS"]
#if idirs:
#    idirs = idirs.split()
#else:
#    idirs = []

#idirs = ['../lib', np.get_include()] #+ idirs

# directories for libraries -L(dir)
#ldirs = os.environ["LDIRS"]
#if ldirs:
#    ldirs = ldirs.split()
#else:
#    ldirs = []

#libs = [] #os.environ['LIBS'].split()

setup(name='spherical_bessel',
      version='0.0.1',
      author='Jun Koda',
      py_modules=['spherical_bessel.sphericial_bessel',
      ],
      ext_modules=[
          Extension('spherical_bessel._spherical_bessel',
                    ['py_package.cpp',
                     'py_spherical_bessel.cpp',
                    ],
                    #include_dirs = idirs,
                    #libraries = libs,
                    undef_macros = ['NDEBUG'],
                    #extra_compile_args = [os.environ['OPT']],
                    #library_dirs = ldirs,

          )
      ],
      packages=['spherical_bessel'],
)
