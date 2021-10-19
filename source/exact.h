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
//   double _recombinations;
// } KwargData;

typedef struct KwargRunResult {
    double r; // not sure what this is
    int recombinations_max;
} KwargRunResult;

extern int exact_randomise;
int beagle(Genes *g, FILE *print_progress);
int beagle_bounded(Genes *g, FILE *print_progress, int lower, int upper);
int beagle_reusable(Genes *g, FILE *print_progress, HashTable *t);
int beagle_reusable_bounded(Genes *g, FILE *print_progress, int lower,
			    int upper, HashTable *t);
LList *beagle_randomised(Genes *g, FILE *print_progress, int r, HashTable *t);
HashTable *beagle_allocate_hashtable(Genes *g, int table_size);
void beagle_deallocate_hashtable(HashTable *t);
double scoring_function(Genes *g);
double score_renormalise(Genes *g, double sc);
KwargRunResult ggreedy(Genes *g, FILE *print_progress, int (*select)(double),
               void (*reset)(void), int ontheflyselection,
               double se_cost, double rm_cost, double r_cost, double rr_cost,
               EList *lookup, int recombinations_max);
#endif
