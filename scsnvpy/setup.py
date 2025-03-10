#Copyright (c) 2018-2020 Gavin W. Wilson
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

from setuptools import setup
from Cython.Build import cythonize
import numpy, os, glob
from os.path import join
from setuptools.extension import Extension
include_dirs = [numpy.get_include()]

extensions = [
    Extension("scsnvpy.data", ["scsnvpy/data.pyx"],
        include_dirs = include_dirs,
    ),
    Extension("scsnvpy.gmix", ["scsnvpy/gmix.pyx"],
        include_dirs = include_dirs,
    ),
    Extension("scsnvpy.snvmats", ["scsnvpy/snvmats.pyx"],
        include_dirs = include_dirs,
        language="c++",
        libraries=['pthread'],
        extra_compile_args=["-std=c++11", "-fPIC"],
        extra_link_args=["-std=c++11"],
    ),
]

pxd_dirs=['scsnvpy/']
setup(name='scsnvpy',
      version='1.0',
      description='scSNV Python Library and Scripts',
      author='Gavin Wilson',
      url='https://github.com/GWW/scsnv',
      author_email='gavin.w.wilson@gmail.com',
      license='MIT',
      packages=['scsnvpy'],
      install_requires=['matplotlib','numpy','scipy', 'cython', 'h5py', 'flammkuchen', 'NCLS', 'tables', 'cyvcf2', 'pandas', 'statsmodels', 'anndata'],
      scripts=['bin/scsnvmisc', 'bin/scsnv2mtx'],
      ext_modules=cythonize(extensions, include_path=pxd_dirs),
)


