#!/usr/bin/env bash

# Run this script in the main folder of SplineLib to update the docker build system

docker build -t mfcats/splinelib:spack -f scripts/docker/Dockerfile.spack .
docker push mfcats/splinelib:spack

docker build -t mfcats/splinelib:gcc8 -f scripts/docker/Dockerfile.gcc8 .
docker push mfcats/splinelib:gcc8