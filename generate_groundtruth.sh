#!/bin/bash

set -e

if ! [ -e build/example/test ]; then
    ./waf configure
    ./waf
fi

mkdir -p ground_truth
cd ground_truth
../build/example/test ../example/test16k.wav

echo "Finished"
