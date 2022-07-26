#!/bin/bash

../arg_builder_source/arg_builder -l -V1 -c -L0 -e1 -r1 -S978175489 -dthreaded_network.dot -R5.0 -B1.3 -M1 < sampleC100.txt

./dot_to_png.sh