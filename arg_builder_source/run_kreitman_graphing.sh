#!/bin/bash

for r in 1.1
do
    for b in 0.5
    do
        ./arg_builder -l -V1 -c -r1 -e2 -ytestyaml.yml -ddotty.dot -L0 -R${r} -B${b} < kreitman_rooted.txt
    done
done