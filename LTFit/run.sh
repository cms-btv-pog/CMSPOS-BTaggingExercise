#!/bin/env sh

WDIR=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${WDIR}/../CFIT

./test
