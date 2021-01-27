#!/usr/bin/env python

"""
This is a setup script for pyJASPAR: a tool for DNA sequence background generation 

This code is free software; you can redistribute it and/or modify it under the terms of the 
BSD License (see the file LICENSE included with the distribution).

@author: Aziz Khan
@email: azez.khan@gmail.com
"""
import os
from distutils.core import setup
from setuptools import find_packages
from pyjaspar import __version__ as VERSION


CLASSIFIERS = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Software Development :: Libraries :: Python Modules',
]

install_requires = [
    'biopython',
]

#def readme():
#    with open('README.rst') as f:
#        return f.read()

def readme(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="pyjaspar",
    description="A serverless interface to Biopython to access different versions of JASPAR database",
    version=VERSION,
    author="Aziz Khan",
    license='GPL',
    platforms='linux/unix',
    author_email="azez.khan@gmail.com",
    url="https://github.com/asntech/pyjaspar",
    long_description=readme("README.rst"),
    long_description_content_type='text/x-rst',
    package_dir={'pyjaspar': 'pyjaspar'},

    packages=['pyjaspar',
        'pyjaspar.data'
        ],

    package_data={'pyjaspar': ['pyjaspar/data/*.sqlite',]},
    include_package_data=True,
    install_requires = install_requires,
    classifiers=CLASSIFIERS,
)
