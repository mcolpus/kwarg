#ifndef _EXACT_H
#define _EXACT_H

#include <stdio.h>
#include "gene.h"
#include "hashtable.h"

// typedef struct _KwargData {
//   double se_cost;
//   double rm_cost;
//   double r_cost;
//   double rr_cost;
//   EList *lookup;
//   int recombinations_max;
//   int rm_max;
//   double g_step_cost;
// } KwargData;

typedef struct KwargRunResult {
    double r; // not sure what this is
    int recombinations_max;
} KwargRunResult;

extern int g_exact_randomise;
int beagle(Genes *g, FILE *print_progress);
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper);
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t);
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
			    int upper, HashTable *t);
LList *beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t);
HashTable *beagle_allocate_hashtable(Genes *g, int table_size);
void beagle_deallocate_hashtable(HashTable *t);
double scoring_function(HistoryFragment *fragment, double temp, double step_cost);
double score_renormalise(HistoryFragment *fragment, double score, double temp, double step_cost);
KwargRunResult ggreedy(PartialHistory *history, FILE *print_progress, int (*select)(double),
               void (*reset)(void),
               double se_cost, double rm_cost, double r_cost, double rr_cost, double temp,
               EList *lookup, int recombinations_max, int print_reference);
#endif
