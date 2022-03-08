#!/bin/bash

python3 ../python_scripts/calculate_optimal_arg_grid.py \
    -i kwarg_records.csv,threaded_records.csv,kwarg_rooted_records.csv,threaded_rooted_records.csv \
    -o kreitman_comparison \
    -l kwarg,threading,rooted_kwarg,rooted_threading