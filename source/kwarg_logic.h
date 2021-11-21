#ifndef KWARG_LOGIC_H
#define KWARG_LOGIC_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

double run_kwarg(Genes *g, FILE *print_progress, int (*select)(double),
               void (*reset)(void), int ontheflyselection, RunSettings run_settings, RunData &run_data);
#endif
