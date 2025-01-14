#ifndef KWARG_LOGIC_H
#define KWARG_LOGIC_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

typedef struct _Result
{
    int seflips = 0;
    int rmflips = 0;
    int recombs = 0;
    int depth;
} Result;

double run_kwarg(Genes *g, FILE *print_progress, int (*select)(double),
                 void (*reset)(void), RunSettings run_settings, RunData &run_data);
std::vector<Result> mass_run_kwarg(Genes *g, FILE *print_progress, std::vector<int> (*take_sample)(std::vector<double>, int),
                                   RunSettings run_settings, RunData &main_path_run_data, int num_samples);
#endif
