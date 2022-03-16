#!/bin/bash

for i in A B C
do
python3 ../python_scripts/graph_results.py -i thread_records_sample${i}25.csv,thread_records_sample${i}50.csv,thread_records_sample${i}100.csv \
                                           -o thread_set_comparison${i}.png \
                                           -l 25,50,100
done