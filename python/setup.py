#!/usr/bin/env python

import os
import sys

if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()

try:
    from setuptools import setup, Extension
    setup, Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension
    setup, Extension

import numpy

include_dirs = [
    numpy.get_include(), "../header", "hodlr",

    # Standard Eigen includes. Provide custom paths with the
    # `--include-dirs=/path/to/eigen3` flag when running `build_ext`.
    "/usr/local/include/eigen3",
    "/usr/local/homebrew/include/eigen3",
    "/opt/local/var/macports/software/eigen3",
    "/opt/local/include/eigen3",
    "/usr/include/eigen3",
]

library_dirs = []
libraries = ["m"]
sources = ["hodlr/_hodlr.cpp", ]
extra_compile_args = ["-ffast-math"]

ext = Extension("hodlr._hodlr",
                sources=sources,
                library_dirs=library_dirs,
                libraries=libraries,
                runtime_library_dirs=library_dirs,
                include_dirs=include_dirs,
                extra_compile_args=extra_compile_args)

setup(
    name="hodlr",
    version="0.3.14",
    author="Daniel Foreman-Mackey",
    author_email="danfm@nyu.edu",
    url="https://github.com/sivaramambikasaran/HODLR_Solver",
    packages=["hodlr"],
    ext_modules=[ext],
    description="HODLR Solver",
    long_description=open("README.rst").read(),
    install_requires=["numpy"],
    classifiers=[
        # "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
