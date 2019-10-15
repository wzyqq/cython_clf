from distutils.core import setup, Extension
from Cython.Build import cythonize

ext = Extension(name="wrap",
				library_dirs=['/home/wzy/gsl/lib'],
                libraries=['m', 'gsl', 'gslcblas'], 
                sources=["clf_mcmc.c", "vars.c","wrap.pyx"])

setup(ext_modules=cythonize(ext))