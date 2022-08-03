#!/bin/bash

for i in 20 40 60 80 100
do
    python3 ../python_scripts/split_into_samples.py -i sampleCall.txt -o sample${i}_ -s ${i}
done
