#! /bin/bash
. /code/SplineLib/scripts/docker/bashrc
spack repo add /code/SplineLib/scripts/spack-repo
cd /code/SplineLib
spack setup splinelib@github
mkdir build
cd build
./../spconfig.py ..
make
/code/SplineLib/build/test/SplineLibTests
