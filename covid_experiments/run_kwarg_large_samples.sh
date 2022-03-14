#!/bin/bash

for j in 1000 all
do
    for i in A B C 
    do
        ../source/kwarg -T10 -Q5 -k -Okwarg_records_sample${i}${j}.csv -S-1,1.5,1.2,1.0,0.8,0.7,0.5,1 \
                                                                       -M-1,1.51,1.21,1.01,0.81,0.71,0.51,1.1 \
                                                                       -R1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,-1 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done