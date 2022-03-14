#!/bin/bash

for j in 50
do
    for i in B 
    do
        ../source/kwarg -T50 -Q20 -k -Otricky_sample${i}${j}.csv -S0.7 -M0.71 -R1 -C2 -r4 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done