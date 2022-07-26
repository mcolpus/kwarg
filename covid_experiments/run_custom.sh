#!/bin/bash

for j in 200
do
    for i in C 
    do
        python3 ../python_scripts/graph_results.py -i thread_records_sample${i}${j}.csv -o thread_sample${i}${j}_graph.png
    done
done
