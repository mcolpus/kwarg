#!/bin/bash

# array=(20)
# #array2=(245 122 81 60 48)
# array2=(245 5 5 5 5)

# for i in "${!array[@]}"; do
#     python3 ../python_scripts/graph_paper_results.py -i records_sample${array[i]} -o graph${array[i]}.png -n ${array2[i]} -c
# done

python3 ../python_scripts/graph_paper_results.py -i records_sample20,records_sample40,records_sample60,records_sample80,records_sample100 -o graph_all.png -n 245,122,81,60,48 -c