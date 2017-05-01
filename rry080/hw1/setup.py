""" setup.py
This code is adapted for 80 columns
|------------------------------------------------------------------------------|
Compiles the cython code into and extension module that can be imported in a
regular python script.
"""
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("ha1utils.pyx",
                          compiler_directives={'embedsignature': True,
                                               'boundscheck': False,
                                               'wraparound': False,
                                               'initializedcheck': False})
)
