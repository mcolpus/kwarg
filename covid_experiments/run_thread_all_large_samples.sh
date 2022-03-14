#!/bin/bash

for j in all
do
    for i in C B
    do
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q5 -R1.9,2.0,3.0,50.0 -B1.3,1.5,2.0,3.0,50.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q5 -R2.0,3.0,50.0 -B1.5,2.0,3.0,50.0 < sample${i}${j}.txt
        echo "completed ${i}${j}"
    done
done
