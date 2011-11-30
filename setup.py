from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "nim",
    ext_modules = cythonize("*.pyx"),
)
