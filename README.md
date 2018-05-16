# SplineLib
Library for spline manipulation.

## License
Copyright 2018 Chair for Computational Analysis of Technical Systems, RWTH Aachen University

SplineLib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation version 3 of the License.

SplineLib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License along with SplineLib.  If not, see
http://www.gnu.org/licenses/.

## Status
[![CodeFactor](https://www.codefactor.io/repository/github/mfcats/splinelib/badge)](https://www.codefactor.io/repository/github/mfcats/splinelib)
[![codecov](https://codecov.io/gh/mfcats/SplineLib/branch/master/graph/badge.svg)](https://codecov.io/gh/mfcats/SplineLib)
[![Build Status](https://travis-ci.org/mfcats/SplineLib.svg?branch=master)](https://travis-ci.org/mfcats/SplineLib)

## Installation

Requirements : fortran compiler, C++11
Advise : Use gcc-Compiler

Clone the repository to your local machine

	 $ git clone https://github.com/mfcats/SplineLib.git

Install spack, see https://spack.readthedocs.io/en/latest/

Add the splinelib repository

	 $ spack repo add path/to/SplineLib/scripts/spack-repo

Now you can install the splinelib

	 $ spack install splinelib

Spack does automatically install the dependencies.

## Installation for developers
Instead of using the install command use

	 $ spack setup splinelib@github

Now you can go to the SplineLib folder and create a build directory in which you can
build the splinelib with cmake and make

	 $ mkdir build
	 $ cd build
	 $ ./setup.py

Testing if required

	 $ ./test/SplineLibTests
