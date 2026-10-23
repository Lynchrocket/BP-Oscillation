from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np


extensions = [
    Extension(
        "decoders.cython.cython_core.bp_core",
        ["src/decoders/cython/cython_core/bp_core.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "decoders.cython.cython_core.osd_core",
        ["src/decoders/cython/cython_core/osd_core.pyx"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        "decoders.cython.cython_core.tesseract_core",
        ["src/decoders/cython/cython_core/tesseract_core.pyx"],
        include_dirs=[np.get_include()],
    ),
]


setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={"language_level": "3"},
    ),
)