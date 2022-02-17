#!/bin/bash

for r in 0.4 0.6 0.9 1.1 1.5 1.9 2.1
do
    for b in 1.1 1.5 1.9 2.1
    do
        ./arg_builder -l -V1 -c -oauto_kreitman_records.csv -Q10 -R${r} -B${b} < kreitman_rooted_ordered.txt
    done
done