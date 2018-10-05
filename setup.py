import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

version = 'x.y.z'
if os.path.exists('VERSION'):
  version = open('VERSION').read().strip()

USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.c'
extensions = [Extension("file_manipulation", ["file_manipulation"+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)
	
setup(
    name='albatradis',
    version=version,
    description='TraDIStron',
	long_description=read('README.md'),
    packages = find_packages(),
    author='Andrew J. Page',
    author_email='andrew.page@quadram.ac.uk',
    url='xxxx',
	ext_modules = extensions,
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
           'biopython >= 1.68',
		   'pyfastaq >= 3.12.0',
		   'scipy',
		   'numpy',
		   'dendropy',
		   'seaborn',
		   'pandas',
		   'graphviz',
		   'cython'
       ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
