#!/bin/bash

# for j in 25 50 100 200
# do
#     for i in A B C 
#     do
#         ../arg_builder_source/arg_builder -l -V0 -c -onew_thread_records_sample${i}${j}.csv -L0 -Q10 -R0.5,0.6,0.9,1.0,1.1,1.5,1.9,2.0,3.0,50.0 -B1.0,1.3,1.5,2.0,3.0,50.0 < sample${i}${j}.txt
#         ../arg_builder_source/arg_builder -l -V0 -c -onew_thread_records_sample${i}${j}.csv -L1 -Q15 -R0.5,0.6,0.9,1.0,1.1,1.5,1.9,2.0,3.0,50.0 -B1.0,1.3,1.5,2.0,3.0,50.0 < sample${i}${j}.txt
#         echo "completed ${i}${j}"
#     done
# done

# for j in 500 1000 all
# do
#     for i in A C B
#     do
#         ../arg_builder_source/arg_builder -l -V0 -c -onew_thread_records_sample${i}${j}.csv -L0 -Q5 -R1.9,2.0,3.0,50.0 -B1.3,1.5,2.0,3.0,50.0 < sample${i}${j}.txt
#         ../arg_builder_source/arg_builder -l -V0 -c -onew_thread_records_sample${i}${j}.csv -L1 -Q5 -R2.0,3.0,50.0 -B1.5,2.0,3.0,50.0 < sample${i}${j}.txt
#         echo "completed ${i}${j}"
#     done
# done

for j in 50
do
    for i in C
    do
        ../source/kwarg -T50,30 -Q10 -k -Orun_kwarg_test_sample${i}${j}.csv -S-1,1.5,1.2,1.1,1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01,1 -M-1,1.51,1.21,1.11,1.01,0.91,0.81,0.71,0.61,0.51,0.41,0.31,0.21,0.11,0.02,1.1 -R1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1 -C2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,-1 < sample${i}${j}.txt
        echo "completed kwarg ${i}${j}"
    done
done