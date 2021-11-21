#ifndef KWARG_LOGIC_H
#define KWARG_LOGIC_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

double scoring_function(Genes *g, double step_cost, RunSettings &run_settings, RunData &run_data);
double score_renormalise(Genes *g, double sc, double step_cost, RunSettings &run_settings, RunData &run_data);
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double),
               void (*reset)(void), int ontheflyselection, RunSettings run_settings, RunData &run_data);
#endif
