#!/bin/bash

for j in 25 50 100 200
do
    for i in A B C 
    do
        python3 ../python_scripts/calculate_optimal_arg_grid.py -i thread_records_sample${i}${j}.csv -o thread_sample${i}${j}
    done
done