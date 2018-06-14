#! /usr/bin/env zsh

# Check command line arguments
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 DIRECTORY" >&2
  exit 1
fi
if ! [ -e "$1" ]; then
  echo "$1 not found" >&2
  exit 1
fi
if ! [ -d "$1" ]; then
  echo "$1 not a directory" >&2
  exit 1
fi

# Get build prefix
BUILD_PREFIX=$1

# Create build folders
mkdir $BUILD_PREFIX/cmake-build-debug-gcc
mkdir $BUILD_PREFIX/cmake-build-release-gcc
mkdir $BUILD_PREFIX/cmake-build-debug-clang
mkdir $BUILD_PREFIX/cmake-build-release-clang
