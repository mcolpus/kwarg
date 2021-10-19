/* common.h: header file for common.c; definition of common data types.
 */
#ifndef _COMMON_H
#define _COMMON_H

#include "llist.h"
#include "gene.h"
#include "arg.h"
/* Function prototypes; a brief explanation of each function should be
 * provided prior to its implementation in common.c.
 */
#ifdef ENABLE_VERBOSE
int verbose();
void set_verbose(int v);
#endif

#ifdef HAPLOTYPE_BLOCKS
typedef struct _SuperColumn {
  int left;
  int right;
} SuperColumn;

extern LList *representativeness;
extern LListCounter *representativeness_counter;
extern int **haploblocks;
void explode_local(int **local, LList *r, int n);
#endif
extern LList *g_eventlist;
extern EList *g_elements;
extern EList *g_sites;
extern int g_seq_numbering;
extern int g_howverbose;
extern double g_step_cost;
extern int g_gene_conversions_enabled;
extern double g_x2random_seed;
extern long int g_xrandom_seed;
extern HashTable *_greedy_functioncalls, *_greedy_beaglereusable;
#ifdef DEBUG
extern HashTable *ancestral_state_trace;
#endif

void *xmalloc(int n);
void *xcalloc(int m, int n);
void *xrealloc(void *oldadr, int n);
#define XRAND_MAX RAND_MAX
void initialise_xrandom();
void initialise_x2random(double seed);
long int xrandom();
long int x2random();
char *i2a(int n);
void pretty_print(FILE *fp, char *s, int l, int i);
void print_option(FILE *fp, char *option, char *description, int l, int i);
void remove_element(int *array, int index, int array_length);
void delete_i(int *array, int i, int array_length);
void delete_by_value(int *array, int v, int array_length);
int max_value(int *array, int array_length);
void print_elist(EList *e, char *comment);
void set_array(double *a1, double *a2, int a2_length, int b);

#endif
