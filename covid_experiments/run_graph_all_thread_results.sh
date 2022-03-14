#!/bin/bash

for j in 25 50 100
do
    for i in A B C 
    do
        python3 ../python_scripts/graph_results.py -i new_thread_records_sample${i}${j}.csv -o new_thread_sample${i}${j}_graph.png -g new_thread_sample${i}${j}_grid.png 
        python3 ../python_scripts/graph_results.py -i new_kwarg_records_sample${i}${j}.csv -o new_kwarg_sample${i}${j}_graph.png

        python3 ../python_scripts/graph_results.py -i new_thread_records_sample${i}${j}.csv,new_kwarg_records_sample${i}${j}.csv \
            -o new_comparison_sample${i}${j}.png -l threaded,kwarg
    done
done