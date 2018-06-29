#! /bin/bash
. /code/bashrc
spack load gcc
spack compiler find
spack install googletest+gmock
