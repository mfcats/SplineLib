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

A compiler with C++11 support is required to build the SplineLib.
The SplineLib was tested with gcc version 7.3.0.

First you have to clone the github repository to your local machine.

	 $ git clone https://github.com/mfcats/SplineLib.git

Spack is used to install all the needed dependencies.
For more information on how to install spack see: https://spack.readthedocs.io/en/latest/

After installing spack you need to add the SplineLib repository.

	 $ spack repo add path/to/SplineLib/scripts/spack-repo

After that you can install the SplineLib.

	 $ spack install splinelib


## Installation for developers

If you wish to install the SplineLib for development use the following command:

	 $ spack setup splinelib@github

Now go to the SplineLib directory and create a build directory in which you can
build the splinelib with cmake and make.

	 $ mkdir build
	 $ cd build
	 $ ./setup.py

To verify that everything went right you can execute the built-in SplineLib tests.

	 $ ./test/SplineLibTests
