/*******************************************************************
 *
 *    tree_kwarg_search.c
 *
 *    Implementation of functions for a heuristic method to
 *    obtain a (near-)minimal ARG in the presence of both recombination
 *    and recurrent mutation (KwARG) using a tree search approach.
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

/* Main function of kwarg implementing neighbourhood search without use of global variables
 */
KwargRunResult thread_safe_kwarg_search(Genes *genes, FILE *print_progress, int (*select_function)(double), void (*reset_select_function)(void),
               double se_cost, double rm_cost, double r_cost, double rr_cost, double temp,
               EList *lookup, int recombinations_max, int print_reference,
               EList *elements, EList *sites)
{
    int i, neighbourhood_size = 0, total_neighbourhood_size = 0, seflips = 0, rmflips = 0, recombs = 0, preds, is_bad_soln = 0;
    double total_event_cost = 0;
    Index *start, *end;
    LList *event_list = MakeLList();
    void (*action)(Genes *);
    const char *names[5];
    names[0] = "Coalescence";
    names[1] = "Sequencing error";
    names[2] = "Recurrent mutation";
    names[3] = "Single recombination";
    names[4] = "Double recombination";
    double *score_array;

    /* Create working copy of genes so as to not change the inputed ones*/
    genes = copy_genes(genes);

    TreeSearchNode *search_node;
    search_node->genes = genes;
    search_node->event_list = event_list;
    search_node->elements = elements;
    search_node->sites = sites;
    search_node->weight = 1; //only doing one run


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
    implode_genes_from(genes, event_list, elements, sites);
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
    /* Store HistoryFragments of possible predecessors in predecessors */
    EList *predecessors;
    predecessors = elist_make();

    int choice_is_fixed = no_recombinations_required(genes);
    if (choice_is_fixed)
        /* Data set can be explained without recombinations */
        free_genes(genes);

    while (!choice_is_fixed)
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
        ac = COAL;
        preds = 0;
        neighbourhood_size = 0;

        _coalesce_compatible_and_entangled_map(genes, action);
        preds = elist_length(_predecessors) - neighbourhood_size;
        neighbourhood_size = elist_length(_predecessors);
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Coalescing entangled: ", preds);
        }

        if (se_cost != -1)
        {
            g_step_cost = se_cost;
            ac = SE;

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
            ac = RM;

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
            ac = RECOMB1;

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
            ac = RECOMB2;

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

        if (!choice_is_fixed)
        {
            score_min = DBL_MAX, score_max = 0;
            for (i = 0; i < elist_length(_predecessors); i++)
            {
                predecessor = (HistoryFragment *)elist_get(_predecessors, i);
                _reset_builtins(predecessor->g); // set predecessor to be _greedy_currentstate
                g_step_cost = predecessor->step_cost;
                // Calculate all the scores and update the min and max
                score_array[i] = scoring_function(predecessor->g, temp, predecessor->step_cost);
            }
        }

        // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
        for (i = 0; i < elist_length(_predecessors); i++)
        {
            predecessor = (HistoryFragment *)elist_get(_predecessors, i);
            _reset_builtins(predecessor->g); // set _greedy_currentstate to be predecessor->g
            
            g_step_cost = predecessor->step_cost;
            double printscore = score_renormalise(predecessor->g, score_array[i], temp, predecessor->step_cost);
            if (print_progress != NULL && g_howverbose == 2)
            {
                fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i + 1, predecessor->step_cost);
                output_genes(predecessor->g, print_progress, NULL);
                print_elist(predecessor->elements, "Sequences: ");
                print_elist(predecessor->sites, "Sites: ");
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
                    free_genes(greedy_choice->g);
                    if (greedy_choice->event != NULL)
                    {
                        while (Length(greedy_choice->event) != 0)
                            free(Pop(greedy_choice->event));
                        DestroyLList(greedy_choice->event);
                    }
                    if (greedy_choice->elements != NULL)
                        elist_destroy(greedy_choice->elements);
                    if (greedy_choice->sites != NULL)
                        elist_destroy(greedy_choice->sites);
                    free(greedy_choice);
                }
                /* Set predecessor to be new choice */
                greedy_choice = predecessor;
            }
            else
            {
                /* Discard predecessor */
                free_genes(predecessor->g);
                if (predecessor->event != NULL)
                {
                    while (Length(predecessor->event) != 0)
                        free(Pop(predecessor->event));
                    DestroyLList(predecessor->event);
                }
                if (predecessor->elements != NULL)
                {
                    elist_destroy(predecessor->elements);
                }
                if (predecessor->sites != NULL)
                {
                    elist_destroy(predecessor->sites);
                }
                free(predecessor);
            }
        }

        free(score_array);

        g_eventlist = tmp;
        elist_empty(_predecessors, NULL); // this should now be empty

        genes = greedy_choice->g;
        g_elements = greedy_choice->elements;
        g_sites = greedy_choice->sites;

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
            Append(g_eventlist, greedy_choice->event);
        }

        total_event_cost += greedy_choice->step_cost;

        /* Clean up */
        free(greedy_choice);
        free(start);
        free(end);
        if (choice_is_fixed)
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
