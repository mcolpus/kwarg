#!/bin/bash

for i in A B C
do
    for j in 25 50 100 200 
    do
        python3 calculate_optimal_arg_grid.py -i sample${i}${j}_records.csv -o graph_sample${i}${j}
        python3 calculate_optimal_arg_grid.py -i randomselection_sample${i}${j}_records.csv -o graph_randomselection_sample${i}${j}
    done
done
