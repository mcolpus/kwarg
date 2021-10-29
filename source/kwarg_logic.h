#ifndef _KWARG_LOGIC_H
#define _KWARG_LOGIC_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

double scoring_function(Genes *g);
double score_renormalise(Genes *g, double sc);
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double),
               void (*reset)(void), int ontheflyselection);
#endif
