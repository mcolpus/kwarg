/* common.h: header file for common.c; definition of common data types.
 */
#ifndef COMMON_H
#define COMMON_H

#include "llist.h"
#include "gene.h"
#include "arg.h"

#define ENABLE_VERBOSE

// These definitions are passed when kwarg is created using make file
#define HAPLOTYPE_BLOCKS
#define BEAGLE_HAPLOTYPEHEURISTIC
#define BEAGLE_DONOTSTORELEAVES
#define ENUMERATE_DONOTSTORELEAVES
#define ENUMERATE_HAPLOTYPEHEURISTIC

#include <vector>

/* Function prototypes; a brief explanation of each function should be
 * provided prior to its implementation in common.c.
 */
#ifdef ENABLE_VERBOSE
int verbose();
void set_verbose(int v);
#endif

#ifdef HAPLOTYPE_BLOCKS
typedef struct _SuperColumn
{
    int left;
    int right;
} SuperColumn;

extern LList *g_representativeness;
extern LListCounter *g_representativeness_counter;
extern int **g_haploblocks;
void explode_local(int **local, LList *r, int n);
#endif

extern LList *g_eventlist;
extern std::vector<int> g_sequence_labels;
extern std::vector<int> g_site_labels;
extern std::vector<int> g_lookup;
extern int g_seq_numbering;
extern double g_se_cost;
extern double rm_cost;
extern double g_r_cost;
extern double g_rr_cost;
extern int g_howverbose;
extern double g_recombinations;
extern int gc_enabled;
extern double g_Temp;
extern double g_run_seed;
extern int g_rec_max, g_rm_max;
extern int g_seed_counter;
extern int g_run_reference;
extern HashTable *g_greedy_functioncalls, *g_greedy_beaglereusable;
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
void print_int_vector(std::vector<int> vec, char *comment);
void set_array(double *a1, double *a2, int a2_length, int b);

/**
 * swaps the values at i and j in the vector
 */
template <typename T>
inline void vector_swap_elements(std::vector<T> &vec, int i, int j)
{
    T temp = std::move(vec[i]);
    vec[i] = std::move(vec[j]);
    vec[j] = std::move(temp);
}

/**
 * append source to destination by moving elements (source left empty)
 */
template <typename T>
inline void vector_append(std::vector<T> &destination, std::vector<T> &source)
{
    if (destination.empty())
        destination = std::move(source);
    else
    {
        std::move(std::begin(source), std::end(source), std::back_inserter(destination));
        // destination.insert(std::end(destination),
        //                    std::make_move_iterator(std::begin(source)),
        //                    std::make_move_iterator(std::end(source)));
    }

}

#endif
