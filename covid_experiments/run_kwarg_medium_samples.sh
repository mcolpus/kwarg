#!/bin/bash

for j in 200 500
do
    for i in A C B 
    do
        ../source/kwarg -T15 -Q10 -k -r600 -Okwarg_records_sample${i}${j}.csv -S-1,1.5,1.2,1.1,1.0,0.9,0.8,0.7,0.5,1 \
                                                                              -M-1,1.51,1.21,1.11,1.01,0.91,0.81,0.71,0.51,1.1 \
                                                                              -R1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,-1 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done