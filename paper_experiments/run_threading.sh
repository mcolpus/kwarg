#!/bin/bash

array=(20 40 60 80 100)
array2=(245 122 81 60 48)
repeat=10

for i in "${!array[@]}"; do
    printf "Starting size ${array[i]} \n"

    for j in $(seq 0 ${array2[i]})
    do
        ../arg_builder_source/arg_builder -l -V1 -c -orecords_sample${array[i]}_${j}.csv -L0 -Q${repeat} -R0.5,0.9,1.0,1.5,2.0,4.0,50.0 -B0.1,0.5,1.0,2.0,50.0 < sample${array[i]}_${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -orecords_sample${array[i]}_${j}.csv -L1 -Q${repeat} -R0.5,0.75,1.0,1.5,4.0,50.0 -B0.3,0.5,1.0,2.0,50.0 < sample${array[i]}_${j}.txt
        
        # Purely recurrent-mutation
        ../arg_builder_source/arg_builder -l -V1 -c -orecords_sample${array[i]}_${j}.csv -L1 -Q${repeat} -R-1 -B0.1,0.5,1.0 < sample${array[i]}_${j}.txt

        # Mostly recomb
        ../arg_builder_source/arg_builder -l -V1 -c -orecords_sample${array[i]}_${j}.csv -L0 -Q${repeat} -M50 -B10 < sample${array[i]}_${j}.txt
        ../arg_builder_source/arg_builder -l -V1 -c -orecords_sample${array[i]}_${j}.csv -L1 -Q${repeat} -M50 -B10 < sample${array[i]}_${j}.txt
        echo "completed ${array[i]} ${j}"
    done
done

