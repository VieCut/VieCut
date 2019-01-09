#!/bin/sh

set -ex

git submodule update --init --recursive

set +x

set -x

# detect number of cores
if [ -e /proc/cpuinfo ]; then
    CORES=$(grep -c ^processor /proc/cpuinfo)
elif [ "$(uname)" == "Darwin" ]; then
    CORES=$(sysctl -n hw.ncpu)
else
    CORES=4
fi
MAKEOPTS=-j$CORES

if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake .. $@
make $MAKEOPTS
ctest -V
