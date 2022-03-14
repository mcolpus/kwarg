#!/bin/bash

for j in 1000 all
do
    for i in A B C 
    do
        # first attempt with rm_cost around 1
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q10 -R0.9,1.0 -B0.1,0.3,0.5,1.0,2.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q10 -R0.9,1.0 -B0.1,0.3,0.5,1.0,2.0 < sample${i}${j}.txt
        # Now try with rm_cost cheaper
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q20 -R0.15,0.3,0.5 -B0.1,0.3,0.5,1.0,2.0,50.0 < sample${i}${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L1 -Q20 -R0.15,0.3,0.5 -B0.1,0.3,0.5,1.0,2.0,50.0 < sample${i}${j}.txt
        # Now look for pure recombination solution
        ../arg_builder_source/arg_builder -l -V1 -c -othread_records_sample${i}${j}.csv -L0 -Q10 -R100.0 -B100.0 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done
