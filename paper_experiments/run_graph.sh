#!/bin/bash

array=(20 40 60 80 100)
array2=(245 122 81 60 48)

for i in "${!array[@]}"; do
    python3 ../python_scripts/graph_paper_results.py -i records_sample${array[i]} -o graph${array[i]}.png -n ${array2[i]}
done
