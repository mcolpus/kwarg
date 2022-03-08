#!/bin/bash

for j in 25 50 100 200
do
    for i in A B C 
    do
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q10 -R0.5,0.6,0.9,1.0,1.1,1.5,1.9,2.0,3.0,50.0 -B1.0,1.3,1.5,2.0,3.0,50.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q10 -R0.5,0.6,0.9,1.0,1.1,1.5,1.9,2.0,3.0,50.0 -B1.0,1.3,1.5,2.0,3.0,50.0 < sample${i}${j}.txt
    done
done
