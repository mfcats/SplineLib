**This repository is no longer maintained. The work is continued under https://github.com/SplineLib/SplineLib.**

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

### Citation
If you want to use SplineLib in your scientific project please cite [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1320259.svg)](https://doi.org/10.5281/zenodo.1320259).

## Status
[![CodeFactor](https://www.codefactor.io/repository/github/mfcats/splinelib/badge)](https://www.codefactor.io/repository/github/mfcats/splinelib)
[![codecov](https://codecov.io/gh/mfcats/SplineLib/branch/master/graph/badge.svg)](https://codecov.io/gh/mfcats/SplineLib)
[![Build Status](https://travis-ci.org/mfcats/SplineLib.svg?branch=master)](https://travis-ci.org/mfcats/SplineLib)


## Installation

A compiler with C++ 17 support is required to successfully build SplineLib. Compilers that fulfill this requirement are GCC 7, GCC 8, clang 5, and clang 6. The project is tested with these compilers. Therefore, the master branch should always run with them.

To use SplineLib, you have to clone the github repository to your local machine.

	 $ git clone https://github.com/mfcats/SplineLib.git

For installation of all dependencies you can use spack.
For more information on how to install spack see: https://spack.readthedocs.io/en/latest/

After installing spack you need to add the SplineLib repository.

	 $ spack repo add path/to/SplineLib/scripts/spack-repo

After that you can install the SplineLib.

	 $ spack install splinelib

This will install the Google testing framework as dependency and afterwards SplineLib. Then, you can use SplineLib by running 

	$ spack load splinelib




## Installation for developers

If you wish to install the SplineLib for development use the following command:

	 $ spack setup splinelib@github

Now go to the SplineLib directory and create a build directory in which you can
build SplineLib with CMake and make.

	 $ mkdir build
	 $ cd build
	 $ ./spconfig.py ..

To verify that everything went right you can execute the built-in SplineLib tests.

	 $ ./test/SplineLibTests
