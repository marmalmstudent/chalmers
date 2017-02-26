from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("ha5utils.pyx",
                          compiler_directives={'embedsignature': True,
                                               'boundscheck': False,
                                               'wraparound': False,
                                               'initializedcheck': False,
                                               'cdivision': True})
)
