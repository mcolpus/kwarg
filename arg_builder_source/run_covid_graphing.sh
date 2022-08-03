#!/bin/bash

for r in 1.1
do
    for b in 0.5
    do
        ./arg_builder -l -V1 -c -e2 -ycovid.yml -dcovid.dot -L0 -R${r} -B${b} < sampleA500.txt
    done
done

python3 add_clade.py