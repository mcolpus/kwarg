#!/bin/bash


python3 ../python_scripts/graph_results.py -i thread_records_sampleA50.csv,thread_records_sampleA100.csv,thread_records_sampleA200.csv,kwarg_records_sampleA50.csv,kwarg_records_sampleA100.csv,kwarg_records_sampleA200.csv \
                                           -o setA_smooth_comparison${i}.png \
                                           -l thread50,thread100,thread200,kwarg50,kwarg100,kwarg200 -t