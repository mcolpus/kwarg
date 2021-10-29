/*******************************************************************
 *
 *    exact.c
 *
 *    Implementation of functions to compute the exact minimum
 *    number of recombinations required under the infinite sites
 *    assumption for an SNP data set (Beagle), and a heuristic method to
 *    obtain a (near-)minimal ARG in the presence of both recombination
 *    and recurrent mutation (KwARG).
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
#include "exact.h"
#include "hashtable.h"
#include "bitfunctions.h"
#include "backtrack.h"
#include "common.h"


/* Update global quantities with contribution from g and free any
 * events that may be stored in g_eventlist.
 */
static int _choice_fixed;
static double _min_ancestral_material, _max_ancestral_material;
static void _update_min_max_ancestral_counts(Genes *g)
{
    int ancestral_material = count_ancestral_material(g);

    if (ancestral_material < _min_ancestral_material)
        _min_ancestral_material = ancestral_material;
    else if (ancestral_material > _max_ancestral_material)
        _max_ancestral_material = ancestral_material;
}

/* Score computation for each state in the neighbourhood */
static double score_min = DBL_MAX, score_max = 0;
double scoring_function(HistoryFragment *history, double temp, double step_cost)
{
    Genes *g = history->g;
    double score;
    double lower_bound;
    int sign;

    // If we have already reached the end, we still cycle through all the possible
    // choices of last step and select the cheapest.
    // We set the score to -(cost of step) if it resolves the last incompatibility,
    // otherwise set the score to -(very big number). The random_select function will
    // pick the move with the least negative score in this case, as needed.
    if (_choice_fixed)
    {
        sign = (temp < 0) - (temp > 0) - (temp == 0);
        if (no_recombinations_required(history_fragment_to_partial_history(history)))
        {
            score = sign * step_cost;
        }
        else
        {
            score = sign * DBL_MAX;
        }
    }
    // If we have not reached the end, score the move as usual.
    else
    {
        // if (_max_ancestral_material < 75)
        // {
        //     // TODO: this uses beagle which is not playing nice
        //     _noexp_rmin();
        //     lower_bound = _greedy_rmin;
        // }
        //         else if(_am < 150) {
        //             lower_bound = _eagl(g);
        //         }
        if (_max_ancestral_material < 200)
        {
            lower_bound = _hb(g);
        }
        else
        {
            lower_bound = hudson_kaplan_genes(g);
        }

        score = (step_cost + lower_bound) * _max_ancestral_material + _am;
        if (score < score_min)
        {
            score_min = score;
        }
        if (score > score_max)
        {
            score_max = score;
        }
    }

    return score;
}

/* Once scores have been computed, renormalise and apply annealing */
double score_renormalise(HistoryFragment *history, double score, double temp, double step_cost)
{
    int sign;

    if (_choice_fixed)
    {
        sign = (temp < 0) - (temp > 0) - (temp == 0);
        if (no_recombinations_required(history_fragment_to_partial_history(history)))
        {
            score = sign * step_cost;
        }
        else
        {
            score = sign * DBL_MAX;
        }
    }
    else
    {
        if (score_max != score_min)
        {
            if (temp != -1)
            {
                score = exp(temp * (1 - (score - score_min) / (score_max - score_min)));
            }
        }
        else
        {
            score = 1;
        }
    }

    return score;
}

/* Store HistoryFragments of possible predecessors in predecessors. It does not deep copy PartialHistory given! */
static EList *_predecessors = NULL;
static Action _ac;
static void _store(PartialHistory *history)
{
    Genes *g = history->g;
    HistoryFragment *history_fragment;

    /* Wrap configuration and events leading to it in a HistoryFragment */
    history_fragment = (HistoryFragment *)xmalloc(sizeof(HistoryFragment));
    history_fragment->event_list = history->event_list;
    history_fragment->g = g;
    history_fragment->step_cost = g_step_cost; //TODO: decide about this
    history_fragment->sequence_labels = history->sequence_labels;
    history_fragment->site_labels = history->site_labels;
    history_fragment->action = _ac;
    if (history->sequence_labels != NULL && g->n != 0 && g->n != elist_length(history->sequence_labels))
    {
        fprintf(stderr, "Error: number of sequence labels in elements [%d] not equal to current size of dataset [%d]. Event cost: %.1f", elist_length(history->sequence_labels), g->n, g_step_cost);
        exit(0);
    }
    if (history->sequence_labels != NULL && g->length > 0 && g->length != elist_length(history->site_labels))
    {
        fprintf(stderr, "Error: number of site labels in sites not equal to current size of dataset.");
        exit(0);
    }
    if (!_choice_fixed && no_recombinations_required(history))
    {
        /* Found a path to the MRCA - choose it */
        _choice_fixed = 1;
    }

    elist_append(_predecessors, history_fragment);
    _update_min_max_ancestral_counts(g);
}

/* Update the lookup list of SE/RM and recombination numbers
 * This is of length recombinations_max, and keeps track of the maximum number of RM events seen
 * for each given number of recombinations already proposed. For example, if we have
 * seen a solution with 5 recombinations and 10 RMs, and the current solution reaches
 * 5 recombinations and 10 RMs but has not yet resolved all incompatibilities, then
 * this solution will be sub-optimal and can be abandoned.
 */
void update_lookup(EList *_lookup, int index, int bd)
{
    int i, j, k;
    // Let S = number of SE + RM in the solution
    // Let R = number of recombinations in the solution
    // Then lookup[R] = S, lookup[R + 1 : R + 2*S] <= S, lookup[R + 2*S : end] = 0
    k = (_lookup->count - 1 > index + 2 * bd ? index + 2 * bd : _lookup->count - 1);
    elist_change(_lookup, index, (void *)bd);
    for (i = index + 1; i <= k; i++)
    {
        j = (int)elist_get(_lookup, i);
        if (j > bd)
        {
            elist_change(_lookup, i, (void *)bd);
        }
    }
    for (i = k + 1; i < _lookup->count; i++)
    {
        elist_change(_lookup, i, (void *)0);
    }
}

/* Main function of kwarg implementing neighbourhood search.
 */
KwargRunResult ggreedy(PartialHistory *history, FILE *print_progress, int (*select_function)(double), void (*reset_select_function)(void),
                       double se_cost, double rm_cost, double r_cost, double rr_cost, double temp,
                       EList *lookup, int recombinations_max, int print_reference)
{
    Genes *genes = history->g;
    int i, neighbourhood_size = 0, total_neighbourhood_size = 0, seflips = 0, rmflips = 0, recombs = 0, preds, is_bad_soln = 0;
    double total_event_cost = 0;
    Index *start, *end;
    LList *tmp = g_eventlist;
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

    if (g_howverbose > 0)
    {
        fprintf(print_progress, "Input data:\n");
        if (g_howverbose == 2)
        {
            output_genes(genes, print_progress, NULL);
        }
        fprintf(print_progress, "%d sequences with %d sites\n", genes->n, genes->length);
    }

    // Reduce the dataset
    implode_genes(history);
    if (g_howverbose > 0)
    {
        printf("%d sequences with %d sites after reducing\n", genes->n, genes->length);
    }
    if (lookup != NULL)
    {
        // set upper bound on run without any recombinations
        if ((int)elist_get(lookup, 0) == INT_MAX)
            update_lookup(lookup, 0, genes->n * genes->length);
    }

    /* Repeatedly choose an event back in time, until data set has been
     * explained.
     */
    _predecessors = elist_make();

    _choice_fixed = no_recombinations_required(history);
    if (_choice_fixed)
        /* Data set can be explained without recombinations */
        // TODO: free history
        free_genes(genes);

    while (!_choice_fixed)
    {
        /* Reset statistics of reachable configurations */
        _min_ancestral_material = INT_MAX;
        _max_ancestral_material = 0;
        HistoryFragment *greedy_choice = NULL;
        neighbourhood_size = 0;
        preds = 0;

        /* Determine interesting recombination ranges */
        start = maximumsubsumedprefixs(genes);
        end = maximumsubsumedpostfixs(genes);

        action = _store;

        /* We have just imploded genes, but we still need to pursue paths
         * coalescing compatible sequences where neither is subsumed in the
         * other but where the ancestral material is still entangled.
         */
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
            fprintf(print_progress, "Searching possible predecessors:\n");
        }
        g_step_cost = 0;
        _ac = COAL;
        preds = 0;
        neighbourhood_size = 0;

        _coalesce_compatible_and_entangled_map(history, action);
        preds = elist_length(_predecessors) - neighbourhood_size;
        neighbourhood_size = elist_length(_predecessors);
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Coalescing entangled: ", preds);
        }

        if (se_cost != -1)
        {
            g_step_cost = se_cost;
            _ac = SE;

            seqerror_flips(genes, action, se_cost);
            preds = elist_length(_predecessors) - neighbourhood_size;
            neighbourhood_size = elist_length(_predecessors);
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Sequencing errors: ", preds);
            }
        }

        if (rm_cost != -1)
        {
            g_step_cost = rm_cost;
            _ac = RM;

            recmut_flips(genes, action, rm_cost);

            preds = elist_length(_predecessors) - neighbourhood_size;
            neighbourhood_size = elist_length(_predecessors);
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Recurrent mutations: ", preds);
            }
        }

        /* Try all sensible events with one split */
        if (r_cost != -1)
        {
            g_step_cost = r_cost;
            _ac = RECOMB1;

            maximal_prefix_coalesces_map(genes, start, end, action);
            preds = elist_length(_predecessors) - neighbourhood_size;
            neighbourhood_size = elist_length(_predecessors);
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Prefix recombinations: ", preds);
            }

            maximal_postfix_coalesces_map(genes, start, end, action);
            preds = elist_length(_predecessors) - neighbourhood_size;
            neighbourhood_size = elist_length(_predecessors);
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Postfix recombinations: ", preds);
            }
        }

        /* Try all sensible events with two splits */
        if (rr_cost != -1)
        {
            g_step_cost = rr_cost;
            _ac = RECOMB2;

            maximal_infix_coalesces_map(genes, start, end, action);
            preds = elist_length(_predecessors) - neighbourhood_size;
            neighbourhood_size = elist_length(_predecessors);
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (infix): ", preds);
            }

            maximal_overlap_coalesces_map(genes, start, end, action);
            preds = elist_length(_predecessors) - neighbourhood_size;
            neighbourhood_size = elist_length(_predecessors);
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (overlap): ", preds);
            }
        }

        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Finished constructing predecessors.", elist_length(_predecessors));
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
        }

        /* Finalise choice and prepare for next iteration */
        free_genes(genes);
        /* Still looking for path to MRCA */
        /* So far we have only enumerated putative predecessors -
         * score these and choose one.
         */

        // Set the tracking lists to NULL for the score computation, and destroy the old elements/sites
        g_eventlist = NULL;
        elist_destroy(g_elements);
        g_elements = NULL;
        elist_destroy(g_sites);
        g_sites = NULL;
        reset_select_function();

        neighbourhood_size = elist_length(_predecessors); // number of predecessors we score
        if (neighbourhood_size == 0)
        {
            fprintf(stderr, "No neighbours left to search but MRCA not reached.");
        }
        total_neighbourhood_size += neighbourhood_size;

        // Calculate all the scores and store in an array
        // Update score_min and score_max for renormalising the score later
        score_array = malloc(elist_length(_predecessors) * sizeof(double));
        HistoryFragment *predecessor;

        if (!_choice_fixed)
        {
            score_min = DBL_MAX, score_max = 0;
            for (i = 0; i < elist_length(_predecessors); i++)
            {
                predecessor = (HistoryFragment *)elist_get(_predecessors, i);
                _reset_builtins(predecessor->g); // set predecessor to be _greedy_currentstate
                g_step_cost = predecessor->step_cost;
                // Calculate all the scores and update the min and max
                score_array[i] = scoring_function(predecessor, 
                                            temp, predecessor->step_cost);
            }
        }

        // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
        for (i = 0; i < elist_length(_predecessors); i++)
        {
            predecessor = (HistoryFragment *)elist_get(_predecessors, i);
            _reset_builtins(predecessor->g); // set _greedy_currentstate to be predecessor->g

            g_step_cost = predecessor->step_cost;
            double printscore = score_renormalise(predecessor, score_array[i], temp, predecessor->step_cost);
            if (print_progress != NULL && g_howverbose == 2)
            {
                fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i + 1, predecessor->step_cost);
                output_genes(predecessor->g, print_progress, NULL);
                print_elist(predecessor->sequence_labels, "Sequences: ");
                print_elist(predecessor->site_labels, "Sites: ");
                fprintf(print_progress, "Predecessor score: %.0f \n\n",
                        (printscore == -DBL_MAX ? -INFINITY : (printscore == DBL_MAX ? INFINITY : printscore)));
                fflush(print_progress);
            }
            if (select_function(printscore))
            {
                // compute score and check if better than that of greedy_choice
                /* If so, discard old choice */
                if (greedy_choice != NULL)
                {
                    free_history_fragment(greedy_choice);
                }
                /* Set predecessor to be new choice */
                greedy_choice = predecessor;
            }
            else
            {
                /* Discard predecessor */
                free_history_fragment(predecessor);
            }
        }

        free(score_array);

        g_eventlist = tmp;
        elist_empty(_predecessors, NULL); // this should now be empty

        genes = greedy_choice->g;
        g_elements = greedy_choice->sequence_labels;
        g_sites = greedy_choice->site_labels;

        switch (greedy_choice->action)
        {
        case COAL:
            break;
        case SE:
            seflips = seflips + greedy_choice->step_cost / se_cost;
            break;
        case RM:
            rmflips = rmflips + greedy_choice->step_cost / rm_cost;
            break;
        case RECOMB1:
            recombs++;
            break;
        case RECOMB2:
            recombs += 2;
            break;
        }

        if (print_progress != NULL && g_howverbose == 2)
        {
            fprintf(print_progress, "%s completed at cost of %.3f.\n", names[greedy_choice->action], greedy_choice->step_cost);
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
            fprintf(print_progress, "Current data:\n");
            output_genes(greedy_choice->g, print_progress, NULL);
            fflush(print_progress);
        }
        if (print_progress != NULL && g_howverbose == 1)
        {
            fprintf(print_progress, "%s at cost %.3f \n", names[greedy_choice->action], greedy_choice->step_cost);
            fflush(print_progress);
        }

        if (g_eventlist != NULL)
        {
            Append(g_eventlist, greedy_choice->event_list);
        }

        total_event_cost += greedy_choice->step_cost;

        /* Clean up */
        free(greedy_choice);
        free(start);
        free(end);
        if (_choice_fixed)
        {
            free_genes(genes);
        }

        // Can abandon the run if the number of recombinations already exceeds recombinations_max
        if (recombs > recombinations_max)
        {
            is_bad_soln = 1;
            break;
        }

        // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
        // seen in earlier solutions.
        if (recombinations_max != INT_MAX && lookup != NULL)
        {
            if (seflips + rmflips > (int)elist_get(lookup, recombs))
            {
                is_bad_soln = 1;
                break;
            }
        }
    }

    // If we exited the loop because of a sub-optimal solution, record this
    if (is_bad_soln)
    {
        fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", print_reference, g_x2random_seed, temp, se_cost, rm_cost, r_cost, rr_cost, total_neighbourhood_size);
    }
    else
    {
        // Otherwise, record the result
        if (print_progress != NULL && g_howverbose > 0)
        {
            fprintf(print_progress, "\nTotal number of states considered: %d\n", total_neighbourhood_size);
            fprintf(print_progress, "Total event cost: %.1f\n", total_event_cost);
            fprintf(print_progress, "%10s %13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Ref", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
                    "SE", "RM", "R", "N_states", "Time");
        }
        fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", print_reference, g_x2random_seed, temp, se_cost, rm_cost, r_cost, rr_cost, seflips, rmflips, recombs, total_neighbourhood_size);

        if (lookup != NULL)
        {
            if (seflips + rmflips < (int)elist_get(lookup, recombs))
            {
                // If found a better bound r < recombinations_max for Rmin, update.
                if (seflips + rmflips == 0)
                {
                    recombinations_max = recombs;
                }
                update_lookup(lookup, recombs, seflips + rmflips);
            }
        }
    }

    elist_destroy(_predecessors);

    KwargRunResult result = {.r = total_event_cost, .recombinations_max = recombinations_max};
    return result;
}
