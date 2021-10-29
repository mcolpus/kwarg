/*******************************************************************
 *
 *    kwarg_logic.c
 *
 *    Implementation of ...
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "gene.h"
#include "common.h"
#include "mergesort.h"
#include "bounds.h"
#include "llist.h"
#include "elist.h"
#include "beagle_logic.h"
#include "kwarg_logic.h"
#include "hashtable.h"
#include "bitfunctions.h"
#include "backtrack.h"
#include "common.h"

/* Functions for interfacing with lower bound computations */
static int _greedy_rmin, _greedy_hk;
static Genes *_greedy_currentstate;
/* Initialise parameters for hashing function call information into m bins */
static int *_greedy_initparam(unsigned long m)
{
    int i, *p = (int *)xmalloc(5 * sizeof(int));

    p[0] = m;
    for (i = 1; i < 5; i++)
        p[i] = xrandom() % m;

    return p;
}

/* Compute minimum number of recombinations needed for current state */
static void *_noexp_rmin()
{
    /* Set up hash table for common use for all beagle invocations */
    if (g_greedy_beaglereusable == NULL) {
        g_greedy_beaglereusable = beagle_allocate_hashtable(_greedy_currentstate, -1);
    }
    if (_greedy_rmin < 0) {
        /* We haven't computed r_min for this configuration yet */
        _greedy_rmin = beagle_reusable(_greedy_currentstate, NULL,
                                       g_greedy_beaglereusable);
    }
    
}


/* Compute haplotype lower bound with the heuristic parameters
 * specified by p.
 */
static int _hb(Genes *g)
{
    int i;
    void **a = (void **)xmalloc(4 * sizeof(void *));

    /* Create array specifying this function call */
    for (i = 0; i < 3; i++) {
        a[i + 1] = (void *)INT_MAX;
    }
    a[0] = (void *)_hb;

    /* Check whether this function call has previously been invoked; if
     * not, compute and store value.
     */
    if (!hashtable_lookup(a, g_greedy_functioncalls, (void **)&i)) {
        i = haplotype_bound_genes(g);
        hashtable_insert(a, (void *)i, g_greedy_functioncalls);
    }

    /* Store lower bound in expression */

    return i;
}

/* Compute lower bound from local exact minimum number of
 * recombinations combined using the composite method.
 */
static int _eagl(Genes *g)
{
    void **a = (void **)xmalloc(2 * sizeof(void *));
    int b, **B;
    Sites *s;

    /* Create array for specifying this function call */
    a[0] = (void *)_eagl;
    a[1] = (void *)(int)(10);

    /* Check whether this function call has previously been invoked; if
     * not, compute and store value.
     */
    if (!hashtable_lookup(a, g_greedy_functioncalls, (void **)&b)) {
        s = genes2sites(g);
        B = hudson_kaplan_local(s);
        free_sites(s);
        b = eagl(g, 10, B, NULL, NULL);
        hashtable_insert(a, (void *)b, g_greedy_functioncalls);
    }

    return b;
}



/* Hash the function call information stored in elm using the
 * parameters in p.
 */
static unsigned long _greedy_hash(void **elm, int *p)
{
    int i;
    unsigned long v = 0;

    if (*elm == _hb)
        for (i = 1; i <= 3; i++)
            v += p[i + 1] * (int)elm[i];
    else
        v = p[1] + p[2] * (int)elm[1];

    return v % p[0];
}

/* Compare the two function calls a and b */
static int _greedy_compare(void **a, void **b)
{
    int i;

    if (*a != *b)
        return 0;

    if (*a == _hb) {
        for (i = 1; i <= 3; i++)
            if ((int)a[i] != (int)b[i])
                return 0;
    }
    else
        return (int)a[1] == (int)b[1];

    return 1;
}

/* Set current ancestral state to g, update am, seq and len to reflect
 * g, and remove information from previous ancestral state from cache.
 */
static double _am, _seq, _len;
static void _reset_builtins(Genes *g)
{
    _greedy_currentstate = g;

    _am = ancestral_material(g);
    _seq = g->n;
    _len = g->length;

    _greedy_rmin = _greedy_hk = -1;
    if (g_greedy_functioncalls != NULL) {
        hashtable_cleanout(g_greedy_functioncalls, free, NULL);
    }
    else
        g_greedy_functioncalls = hashtable_new
                                (6, (unsigned long (*)(void *, void *))_greedy_hash,
                                 (int (*)(void *, void *))_greedy_compare,
                                 (void *(*)(unsigned long))_greedy_initparam);
}

/* Update global quantities with contribution from g and free any
 * events that may be stored in g_eventlist.
 */
static int _choice_fixed;
HistoryFragment *_greedy_choice;
static double _minam, _maxam, _minseq, _maxseq, _minlen, _maxlen;
static Action ac;
static void __update(Genes *g)
{
    int am = ancestral_material(g);

    if (am < _minam)
        _minam = am;
    else if (am > _maxam)
        _maxam = am;
    if (g->n < _minseq)
        _minseq = g->n;
    else if (g->n > _maxseq)
        _maxseq = g->n;
    if (g->length < _minlen)
        _minlen = g->length;
    else if (g->length > _maxlen)
        _maxlen = g->length;
}

static void _update(Genes *g)
{
    /* Check whether we have found a path to the MRCA */
    if (!_choice_fixed && no_recombinations_required(g)) {
        /* Found a path to the MRCA - choose it */
        _choice_fixed = 1;
        _greedy_choice = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));
        _greedy_choice->g = g;
        _greedy_choice->event = g_eventlist;
        _greedy_choice->recombinations = g_step_cost;
        _greedy_choice->elements = g_sequence_labels;
        _greedy_choice->sites = g_site_labels;
        _greedy_choice->action = ac;
        
        return;
    }

    /* Update global quantities */
    if (!_choice_fixed)
        __update(g);

    /* Free memory used for g and event leading to it */
    if (g_eventlist != NULL) {
        while (Length(g_eventlist) != 0)
            free(Pop(g_eventlist));
        DestroyLList(g_eventlist);
    }
    free_genes(g);
}

/* Score computation for each state in the neighbourhood */
static double sc_min = DBL_MAX, sc_max = 0;
static double prev_lb = 0, current_lb = 0, _lb;
double scoring_function(Genes *g) {
    
    double sc;
    double lb;
    int sign;
    
    // If we have already reached the end, we still cycle through all the possible
    // choices of last step and select the cheapest.
    // We set the score to -(cost of step) if it resolves the last incompatibility,
    // otherwise set the score to -(very big number). The random_select function will
    // pick the move with the least negative score in this case, as needed.
    if(_choice_fixed) {
        sign = (g_Temp < 0) - (g_Temp > 0) - (g_Temp == 0);
        if(no_recombinations_required(g)) {
            sc = sign * g_step_cost;
        }
        else {
            sc = sign * DBL_MAX;
        }
    }
    // If we have not reached the end, score the move as usual. 
    else {
        if(_maxam < 75) {
            _noexp_rmin();
            lb = _greedy_rmin;
        }
//         else if(_am < 150) {
//             lb = _eagl(g);
//         }
        else if(_maxam < 200){
            lb = _hb(g);
        }
        else {
            lb = hudson_kaplan_genes(g);
        }
        
        _lb = lb;
        
        sc = (g_step_cost + lb) * _maxam + _am;
        if(sc < sc_min) {
            sc_min = sc;
        }
        if(sc > sc_max) {
            sc_max = sc;
        }
    }
    
    return sc;
}

/* Once scores have been computed, renormalise and apply annealing */
double score_renormalise(Genes *g, double sc) {
    
    int sign;
    
    if(_choice_fixed) {
        sign = (g_Temp < 0) - (g_Temp > 0) - (g_Temp == 0);
        if(no_recombinations_required(g)) {
            sc = sign * g_step_cost;
        }
        else {
            sc = sign * DBL_MAX;
        }
    }
    else {
        if(sc_max != sc_min) {
            if(g_Temp != -1) {
                sc = exp(g_Temp * (1 - (sc - sc_min)/(sc_max - sc_min)));
            } 
        }
        else {
            sc = 1;
        }
    }
    
    return sc;
}

/* Store HistoryFragments of possible predecessors in predecessors */
static EList *_predecessors = NULL;
static void _store(Genes *g)
{
    HistoryFragment *f;
    
    /* Wrap configuration and events leading to it in a HistoryFragment */
    f = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));
    f->event = g_eventlist;
    f->g = g;
    f->recombinations = g_step_cost;
    f->elements = g_sequence_labels;
    f->sites = g_site_labels;
    f->action = ac;
    if(g_sequence_labels != NULL && g->n != 0 && g->n != elist_length(g_sequence_labels)) {
        fprintf(stderr, "Error: number of sequence labels in g_sequence_labels [%d] not equal to current size of dataset [%d]. Event type: %.1f", elist_length(g_sequence_labels), g->n, g_step_cost);
        exit(0);
    }
    if(g_sequence_labels != NULL && g->length > 0 && g->length != elist_length(g_site_labels)) {
        fprintf(stderr, "Error: number of site labels in g_site_labels not equal to current size of dataset.");
        exit(0);
    }
    if (!_choice_fixed && no_recombinations_required(g)) {
        /* Found a path to the MRCA - choose it */
        _choice_fixed = 1;
//         _greedy_choice = f;
    }

    elist_append(_predecessors, f);
    __update(g);

}

static int (*_choice_function)(double);


/* Update the lookup list of SE/RM and recombination numbers
 * This is of length g_recombs_max, and keeps track of the maximum number of RM events seen
 * for each given number of recombinations already proposed. For example, if we have
 * seen a solution with 5 recombinations and 10 RMs, and the current solution reaches
 * 5 recombinations and 10 RMs but has not yet resolved all incompatibilities, then
 * this solution will be sub-optimal and can be abandoned.
 */
void update_lookup(EList *lku, int index, int bd) {
    int i, j, k;
    // Let S = number of SE + RM in the solution
    // Let R = number of recombinations in the solution
    // Then lookup[R] = S, lookup[R + 1 : R + 2*S] <= S, lookup[R + 2*S : end] = 0
    k = (lku->count - 1 > index + 2*bd ? index + 2*bd : lku->count - 1);
    elist_change(lku, index, (void *)bd);
    for(i = index + 1; i <= k; i++) {
        j = (int)elist_get(lku, i);
        if(j > bd) {
            elist_change(lku, i, (void *)bd);
        }
    }
    for(i = k+1; i < lku->count; i++) {
        elist_change(lku, i, (void *)0);
    }
}

/* Main function of kwarg implementing neighbourhood search.
 */
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double), void (*reset)(void), int ontheflyselection)
{
    int global, i, nbdsize = 0, total_nbdsize = 0, seflips = 0, rmflips = 0, recombs = 0, preds, bad_soln = 0;
    double r = 0;
    Index *start, *end;
    LList *tmp = g_eventlist;
    double printscore = 0;
    HistoryFragment *f;
    void (*action)(Genes *);
    const char *names[5];
    names[0] = "Coalescence";
    names[1] = "Sequencing error";
    names[2] = "Recurrent mutation";
    names[3] = "Single recombination";
    names[4] = "Double recombination";
    double *score_array;
    
    #ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
    #endif
    
    if(g_rm_max < INT_MAX) {
        update_lookup(g_lookup, 0, g_rm_max);
    }
    
    /* Create working copy of g */
    g = copy_genes(g);
    
    if(g_howverbose > 0) {
        fprintf(print_progress, "Input data:\n");
        if(g_howverbose == 2) {
            output_genes(g, print_progress, NULL);
            print_elist(g_sequence_labels, "g_elements: ");
            print_elist(g_site_labels, "g_sites: ");
        }
        fprintf(print_progress, "%d sequences with %d sites\n", g->n, g->length);
    }

    // Reduce the dataset
    implode_genes(g);
    if(g_howverbose > 0) {
        printf("%d sequences with %d sites after reducing\n", g->n, g->length);
    }
    if(g_lookup != NULL) {
        if((int)elist_get(g_lookup, 0) == INT_MAX)
            update_lookup(g_lookup, 0, g->n * g->length);
    }

    global = 1;
    
    /* Repeatedly choose an event back in time, until data set has been
     * explained.
     */
    _choice_function = select;
    if (!ontheflyselection && global)
        _predecessors = elist_make();
    if ((_choice_fixed = no_recombinations_required(g)) != 0)
        /* Data set can be explained without recombinations */
        free_genes(g);
    
     while (!_choice_fixed) {
        /* Reset statistics of reachable configurations */
        _minam = _minseq = _minlen = INT_MAX;
        _maxam = _maxseq = _maxlen = 0;
        _greedy_choice = NULL;
        nbdsize = 0;
        preds = 0;

        /* Determine interesting recombination ranges */
        start = maximumsubsumedprefixs(g);
        end = maximumsubsumedpostfixs(g);
     

        action = _store;
        
        /* We have just imploded genes, but we still need to pursue paths
            * coalescing compatible sequences where neither is subsumed in the
            * other but where the ancestral material is still entangled.
            */
        if(g_howverbose > 0) {
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
            fprintf(print_progress, "Searching possible predecessors:\n");
        }
        g_step_cost = 0;
        ac = COAL;
        preds = 0;
        nbdsize = 0;
        
        coalesce_compatibleandentangled_map(g, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(g_howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Coalescing entangled: ", preds);
            }
        
        if(g_se_cost != -1) {
            g_step_cost = g_se_cost;
            ac = SE;
            
            seqerror_flips(g, action);
                preds = elist_length(_predecessors) - nbdsize;
                nbdsize = elist_length(_predecessors);
                if(g_howverbose > 0) {
                    fprintf(print_progress, "%-40s %3d\n", "Sequencing errors: ", preds);
                }
        }
        
        if(g_rm_cost != -1) {
            g_step_cost = g_rm_cost;
            ac = RM;
            
            recmut_flips(g, action);

            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(g_howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Recurrent mutations: ", preds);
            }
        }
        
        /* Try all sensible events with one split */
        if(g_r_cost != -1) {
            g_step_cost = g_r_cost;
            ac = RECOMB1;
            
            maximal_prefix_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(g_howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Prefix recombinations: ", preds);
            }
            
            maximal_postfix_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(g_howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Postfix recombinations: ", preds);
            }
        }
            
        /* Try all sensible events with two splits */
        if(g_rr_cost != -1) {
            g_step_cost = g_rr_cost;
            ac = RECOMB2;
            
            maximal_infix_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(g_howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (infix): ", preds);
            }
            
            maximal_overlap_coalesces_map(g, start, end, action);
            preds = elist_length(_predecessors) - nbdsize;
            nbdsize = elist_length(_predecessors);
            if(g_howverbose > 0) {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (overlap): ", preds);
            }
        }
        
        if(g_howverbose > 0) {
            fprintf(print_progress, "%-40s %3d\n", "Finished constructing predecessors.", elist_length(_predecessors));
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
        }
        
        
        /* Finalise choice and prepare for next iteration */
        free_genes(g);
            /* Still looking for path to MRCA */
            if (!ontheflyselection && global) {
                /* So far we have only enumerated putative predecessors -
                 * score these and choose one.
                 */
                
                // Set the tracking lists to NULL for the score computation, and destroy the old elements/sites
                g_eventlist = NULL;
                elist_destroy(g_sequence_labels);
                g_sequence_labels = NULL;
                elist_destroy(g_site_labels);
                g_site_labels = NULL;
                reset();
                
                nbdsize = elist_length(_predecessors); // number of predecessors we score
                if(nbdsize == 0) {
                    fprintf(stderr, "No neighbours left to search but MRCA not reached.");
                }
                total_nbdsize = total_nbdsize + nbdsize;
                
                // Calculate all the scores and store in an array
                // Update sc_min and sc_max for renormalising the score later
                score_array = malloc(elist_length(_predecessors) * sizeof(double));
                if(!_choice_fixed) {
                    sc_min = DBL_MAX, sc_max = 0;
                    for (i = 0; i < elist_length(_predecessors); i++) {
                        f = (HistoryFragment *)elist_get(_predecessors, i);
                        _reset_builtins(f->g); // set f to be _greedy_currentstate
                        g_step_cost = f->recombinations;
                        // Calculate all the scores and update the min and max
                        score_array[i] = scoring_function(f->g);
                    }
                }
                
                // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
                for (i = 0; i < elist_length(_predecessors); i++) {
                    f = (HistoryFragment *)elist_get(_predecessors, i);
                    _reset_builtins(f->g); // set _greedy_currentstate to be f->g
                    // Bug fix: need to update g_step_cost otherwise this will always be 2
                    g_step_cost = f->recombinations; 
                    printscore = score_renormalise(f->g, score_array[i]);
                    if (print_progress != NULL && g_howverbose == 2) {
                        fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i+1, f->recombinations);
                        output_genes(f->g, print_progress, NULL);
                        print_elist(f->elements, "Sequences: ");
                        print_elist(f->sites, "Sites: ");
                        fprintf(print_progress, "Predecessor score: %.0f \n\n", 
                                (printscore == -DBL_MAX ? -INFINITY : (printscore == DBL_MAX ? INFINITY : printscore)));
                        fflush(print_progress);
                    }
                    if (select(printscore)) {
                        // compute score and check if better than that of _greedy_choice
                        /* If so, discard old choice */
                        if (_greedy_choice != NULL) {
                            free_genes(_greedy_choice->g);
                            if (_greedy_choice->event != NULL) {
                                while (Length(_greedy_choice->event) != 0)
                                    free(Pop(_greedy_choice->event));
                                DestroyLList(_greedy_choice->event);
                            }
                            if(_greedy_choice->elements != NULL)
                                elist_destroy(_greedy_choice->elements);
                            if(_greedy_choice->sites != NULL)
                                elist_destroy(_greedy_choice->sites);
                            free(_greedy_choice);
                        }
                        /* Set f to be new choice */
                        _greedy_choice = f;
                    }
                    else {
                            /* Discard f */
                            free_genes(f->g);
                            if (f->event != NULL) {
                                while (Length(f->event) != 0)
                                    free(Pop(f->event));
                                DestroyLList(f->event);
                            }
                            if(f->elements != NULL) {
                                elist_destroy(f->elements);
                            }
                            if(f->sites != NULL) {
                                elist_destroy(f->sites);
                            }
                            free(f);
                    }
                }
                
                free(score_array);
                
                g_eventlist = tmp;
                elist_empty(_predecessors, NULL); // this should now be empty
            }
            
            g = _greedy_choice->g;
            g_sequence_labels = _greedy_choice->elements;
            g_site_labels = _greedy_choice->sites;
            
            switch(_greedy_choice->action) {
                case COAL:
                    break;
                case SE:
                    seflips = seflips + _greedy_choice->recombinations/g_se_cost;
                    break;
                case RM:
                    rmflips = rmflips + _greedy_choice->recombinations/g_rm_cost;
                    break;
                case RECOMB1:
                    recombs++;
                    break;
                case RECOMB2:
                    recombs += 2;
                    break;
            }
            
            if (print_progress != NULL && g_howverbose == 2) {
                fprintf(print_progress, "%s completed at cost of %.3f.\n", names[_greedy_choice->action], _greedy_choice->recombinations);
                fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
                fprintf(print_progress, "Current data:\n");
                output_genes(_greedy_choice->g, print_progress, NULL);
                fflush(print_progress);
            }
            if (print_progress != NULL && g_howverbose == 1) {
                fprintf(print_progress, "%s at cost %.3f \n", names[_greedy_choice->action], _greedy_choice->recombinations);
                fflush(print_progress);
            }
            /* Predecessor and events leading to it are stored in _greedy_choice */
//         }
        
        
        if (g_eventlist != NULL) {
            Append(g_eventlist, _greedy_choice->event);
        }
        
        r += _greedy_choice->recombinations;
        
        /* Clean up */
        free(_greedy_choice);
        free(start);
        free(end);
        if(_choice_fixed) {
            free_genes(g);
        }
        
        // Can abandon the run if the number of recombinations already exceeds g_recombs_max
        if(recombs > g_recombs_max) {
            bad_soln = 1;
            break;
        }
        
        // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
        // seen in earlier solutions.
        if(g_recombs_max != INT_MAX && g_lookup != NULL) {
            if(seflips + rmflips > (int)elist_get(g_lookup, recombs)) {
                bad_soln = 1;
                break;
            }
        }
        
    }
    
    // If we exited the loop because of a sub-optimal solution, record this
    if(bad_soln) {
        if(g_reference > 0){
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", g_reference, g_r_seed, g_Temp, g_se_cost, g_rm_cost, g_r_cost, g_rr_cost, total_nbdsize);
        } else {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", g_r_seed, g_Temp, g_se_cost, g_rm_cost, g_r_cost, g_rr_cost, total_nbdsize);
        }
    }
    else {
    // Otherwise, record the result
        if (print_progress != NULL && g_howverbose > 0) {
            fprintf(print_progress, "\nTotal number of states considered: %d\n", total_nbdsize);
            fprintf(print_progress, "Total event cost: %.1f\n", r);
            if(g_reference > 0) {
                fprintf(print_progress, "%10s %13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Ref", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
                    "SE", "RM", "R", "N_states", "Time");
            } else {
                fprintf(print_progress, "%13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost", "SE", "RM", "R", "N_states", "Time");
            }
        }
        if(g_reference > 0) {
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", g_reference, g_r_seed, g_Temp, g_se_cost, g_rm_cost, g_r_cost, g_rr_cost, seflips, rmflips, recombs, total_nbdsize);
        } else {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", g_r_seed, g_Temp, g_se_cost, g_rm_cost, g_r_cost, g_rr_cost, seflips, rmflips, recombs, total_nbdsize);
        }
        if(g_lookup != NULL) {
            if(seflips + rmflips < (int)elist_get(g_lookup, recombs)){
                // If found a better bound r < g_recombs_max for Rmin, update.
                if(seflips + rmflips == 0) {
                    g_recombs_max = recombs;
                }
                update_lookup(g_lookup, recombs, seflips + rmflips);
            }
        }
    }
    
    elist_destroy(_predecessors);
    

    
    return r;

}

