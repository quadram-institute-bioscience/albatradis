import glob
import os

from setuptools import setup, find_packages, Extension

requirements = [line.rstrip() for line in open("requirements.txt", "rt")]

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
    description='Comparitive transposon mutagenesis experiment analysis',
	long_description='AlbaTraDIS is a software application for performing rapid large-scale comparative analysis of TraDIS experiments (transposon mutagenesis) whilst also predicting the impact of inserts on nearby genes. It allows for experiements with multiple conditions to be easily analysed using statistical methods developed in the QuaTraDIS toolkit.',
    long_description_content_type='text/markdown',
    packages = find_packages(),
    author='Andrew J. Page',
    author_email='andrew.page@quadram.ac.uk',
    url='https://github.com/quadram-institute-bioscience/albatradis',
	ext_modules = extensions,
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=requirements,
    extras_require={
        "dev": ["semantic-version", "pytest-cov"],
        "test": ["pytest-cov"]
    },
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
		'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)'
    ],
)
