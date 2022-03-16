#!/bin/bash

python3 ../python_scripts/split_records.py -i split_test_records.csv -k 25 -g

python3 ../python_scripts/graph_results.py -i split_test_records_1.csv,split_test_records_2.csv,split_test_records_3.csv,split_test_records_4.csv,split_test_records_5.csv \
            -o split_test_comparison.png -l 1,2,3,4,5

for i in {1..25}
do
rm split_test_records_${i}.csv
done
