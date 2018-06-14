#!/usr/bin/env bash

if ! which spack >/dev/null; then
    mkdir -p $SPACK_ROOT
    git clone --depth 50 https://github.com/spack/spack.git $SPACK_ROOT
    echo -e "config:""\n  build_jobs:"" 2" > $SPACK_ROOT/etc/spack/config.yaml;
    spack bootstrap
fi