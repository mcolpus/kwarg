#!/bin/bash

for name in kwarg kwarg_rooted kwargDFS kwargDFS_rooted threaded threaded_rooted
do
    python3 calculate_optimal_arg_grid.py -i ${name}_records.csv -o ${name}
done

python3 calculate_optimal_arg_grid.py -i kwarg_rooted_records.csv,kwargDFS_rooted_records.csv -o kwarg_comparison

python3 calculate_optimal_arg_grid.py -i kwarg_rooted_records.csv,threaded_rooted_records.csv -o threading_comparison
