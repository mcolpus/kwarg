#ifndef EXACT_H
#define EXACT_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

extern int exact_randomise;
int beagle(Genes *g, FILE *print_progress);
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper);
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t);
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
			    int upper, HashTable *t);
EventList beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t);
HashTable *beagle_allocate_hashtable(Genes *g, int table_size);
void beagle_deallocate_hashtable(HashTable *t);
double scoring_function(Genes *g, double step_cost, RunSettings &run_settings);
double score_renormalise(Genes *g, double sc, double step_cost, RunSettings &run_settings);
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double),
               void (*reset)(void), int ontheflyselection, RunSettings run_settings);
#endif
