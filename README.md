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

Requirements : fortran compiler
Prefered : all gnu compiler

git clone https://github.com/mfcats/SplineLib.git
cd SplineLib

install spack (Do not forget about the .zshrc-file)
spack bootstrap
spack repo add path/to/SplineLib/scripts/spack-repo
spack setup splinelib@github -> this requires a fortran compiler
spack load cmake


# usefull commands for spack
spack spec splinelib -> to see dependencies
spack location -i cmake -> to see location

# build project
mkdir build
cd build
$GTEST_ROOT = spack location -i googletest+gmock
cmake -DGTEST_ROOT=$GTEST_ROOT ..
make 
