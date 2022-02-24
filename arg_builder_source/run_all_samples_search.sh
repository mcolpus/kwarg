#!/bin/bash

for i in C
do
    for j in 25 50 100 200 500 
    do
        ./arg_builder -l -V1 -c -osample${i}${j}_records.csv -r1 -R0.4,0.6,0.9,1.0,1.1,1.5,1.9,2.1,2.5,3.0,50.0 -B1.1,1.3,1.5,1.9,2.1,2.9,3.1,50.0 -Q10 < sampleref${i}${j}.txt
    done
done
