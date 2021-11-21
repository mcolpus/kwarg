/*******************************************************************
 *
 *    kwarg_logic.cpp
 *
 *    Implementation of functions to compute a heuristic method to
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

#include <memory>

#include "gene.h"
#include "common.h"
#include "mergesort.h"
#include "bounds.h"
#include "llist.h"
#include "kwarg_logic.h"
#include "beagle_logic.h"
#include "hashtable.h"
#include "bitfunctions.h"
#include "backtrack.h"
#include "common.h"


/* Update global quantities with contribution from g
 */
static int _choice_fixed;
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

/* Score computation for each state in the neighbourhood */
static double sc_min = DBL_MAX, sc_max = 0;
static double prev_lb = 0, current_lb = 0, _lb;
double scoring_function(Genes *g, double step_cost, RunSettings &run_settings, RunData &run_data)
{
    double sc;
    double lb;
    int sign;

    // If we have already reached the end, we still cycle through all the possible
    // choices of last step and select the cheapest.
    // We set the score to -(cost of step) if it resolves the last incompatibility,
    // otherwise set the score to -(very big number). The random_select function will
    // pick the move with the least negative score in this case, as needed.
    if (_choice_fixed)
    {
        sign = (run_settings.temp < 0) - (run_settings.temp > 0) - (run_settings.temp == 0);
        if (no_recombinations_required(g, run_data))
        {
            sc = sign * step_cost;
        }
        else
        {
            sc = sign * DBL_MAX;
        }
    }
    // If we have not reached the end, score the move as usual.
    else
    {
        if (_maxam < 75)
        {
            lb = noexp_rmin(run_data);
        }
        //         else if(_am < 150) {
        //             lb = _eagl(g);
        //         }
        else if (_maxam < 200)
        {
            lb = hb(g, run_data);
        }
        else
        {
            lb = hudson_kaplan_genes(g);
        }

        _lb = lb;

        sc = (step_cost + lb) * _maxam + get_am();
        if (sc < sc_min)
        {
            sc_min = sc;
        }
        if (sc > sc_max)
        {
            sc_max = sc;
        }
    }
    return sc;
}

/* Once scores have been computed, renormalise and apply annealing */
double score_renormalise(Genes *g, double sc, double step_cost, RunSettings &run_settings, RunData &run_data)
{

    int sign;

    if (_choice_fixed)
    {
        sign = (run_settings.temp < 0) - (run_settings.temp > 0) - (run_settings.temp == 0);
        if (no_recombinations_required(g, run_data))
        {
            sc = sign * step_cost;
        }
        else
        {
            sc = sign * DBL_MAX;
        }
    }
    else
    {
        if (sc_max != sc_min)
        {
            if (run_settings.temp != -1)
            {
                sc = exp(run_settings.temp * (1 - (sc - sc_min) / (sc_max - sc_min)));
            }
        }
        else
        {
            sc = 1;
        }
    }

    return sc;
}

static int (*_choice_function)(double);

/* Update the lookup list of SE/RM and recombination numbers
 * This is of length rec_max, and keeps track of the maximum number of RM events seen
 * for each given number of recombinations already proposed. For example, if we have
 * seen a solution with 5 recombinations and 10 RMs, and the current solution reaches
 * 5 recombinations and 10 RMs but has not yet resolved all incompatibilities, then
 * this solution will be sub-optimal and can be abandoned.
 */
void update_lookup(std::vector<int> &lku, int index, int bd)
{
    int i, k;
    // Let S = number of SE + RM in the solution
    // Let R = number of recombinations in the solution
    // Then lookup[R] = S, lookup[R + 1 : R + 2*S] <= S, lookup[R + 2*S : end] = 0
    k = (lku.size() - 1 > index + 2 * bd ? index + 2 * bd : lku.size() - 1);
    lku[index] = bd;
    for (i = index + 1; i <= k; i++)
    {
        if (lku[i] > bd)
        {
            lku[i] = bd;
        }
    }
    for (i = k + 1; i < lku.size(); i++)
    {
        lku[i] = 0;
    }
}

/* Main function of kwarg implementing neighbourhood search.
 */
double ggreedy(Genes *g, FILE *print_progress, int (*select)(double), void (*reset)(void), int ontheflyselection,
               RunSettings run_settings, RunData &main_path_run_data)
{
    int global, nbdsize = 0, total_nbdsize = 0, seflips = 0, rmflips = 0, recombs = 0, preds, bad_soln = 0;
    double r = 0;
    Index *start, *end;
    double printscore = 0;
    const char *names[5];
    names[0] = "Coalescence";
    names[1] = "Sequencing error";
    names[2] = "Recurrent mutation";
    names[3] = "Single recombination";
    names[4] = "Double recombination";
    double *score_array;

    /* Store HistoryFragments of possible predecessors in predecessors */
    auto predecessors = std::vector<std::unique_ptr<HistoryFragment>>{};

#ifdef ENABLE_VERBOSE
    int v = verbose();
    set_verbose(0);
#endif

    if (run_settings.rm_max < INT_MAX)
    {
        update_lookup(g_lookup, 0, run_settings.rm_max);
    }

    /* Create working copy of g */
    g = copy_genes(g);

    if (g_howverbose > 0)
    {
        fprintf(print_progress, "Input data:\n");
        if (g_howverbose == 2)
        {
            output_genes(g, print_progress, NULL);
        }
        fprintf(print_progress, "%d sequences with %d sites\n", g->n, g->length);
    }

    // Reduce the dataset
    implode_genes(g, main_path_run_data);
    if (g_howverbose > 0)
    {
        if (g_howverbose == 2)
        {
            output_genes(g, print_progress, NULL);
        }
        printf("%d sequences with %d sites after reducing\n", g->n, g->length);
    }
    if (!g_lookup.empty())
    {
        if (g_lookup[0] == INT_MAX)
            update_lookup(g_lookup, 0, g->n * g->length);
    }

    global = 1;

    /* Repeatedly choose an event back in time, until data set has been
     * explained.
     */
    _choice_function = select;
    if (!ontheflyselection && global)
        predecessors.clear();
    _choice_fixed = no_recombinations_required(g, main_path_run_data);
    if (_choice_fixed != 0)
        /* Data set can be explained without recombinations */
        free_genes(g);

    // greedy_choice will point to History fragment to be selected
    auto greedy_choice = std::make_unique<HistoryFragment>();

    while (!_choice_fixed)
    {
        /* Reset statistics of reachable configurations */
        _minam = _minseq = _minlen = INT_MAX;
        _maxam = _maxseq = _maxlen = 0;
        // g = copy_genes(g);
        greedy_choice.reset(nullptr); // = NULL;
        nbdsize = 0;
        preds = 0;

        /* Determine interesting recombination ranges */
        start = maximumsubsumedprefixs(g);
        end = maximumsubsumedpostfixs(g);

        // Store current state to a HistoryFragment
        // run_data takes parameter by copy!
        auto action = [&](Genes *g, RunData run_data)
        {
            auto f = std::make_unique<HistoryFragment>();

            /* Wrap configuration and events leading to it in a HistoryFragment */
            f->events = run_data.eventlist;
            f->g = g;
            f->step_cost = run_data.current_step_cost;
            f->elements = std::move(run_data.sequence_labels);
            f->sites = std::move(run_data.site_labels);
            f->action = ac;
            if (!f->elements.empty() && g->n != 0 && g->n != f->elements.size())
            {
                fprintf(stderr, "Error: number of sequence labels in sequence_labels [%d] not equal to current size of dataset [%d]. Event type: %.1f", f->elements.size(), g->n, run_data.current_step_cost);
                exit(0);
            }
            if (!f->elements.empty() && g->length > 0 && g->length != f->sites.size())
            {
                fprintf(stderr, "Error: number of site labels in sites not equal to current size of dataset.");
                exit(0);
            }
            if (!_choice_fixed && no_recombinations_required(g, run_data))
            {
                /* Found a path to the MRCA - choose it */
                _choice_fixed = 1;
                //         greedy_choice = f;
            }

            predecessors.push_back(std::move(f));
            __update(g);
        };

        /* We have just imploded genes, but we still need to pursue paths
         * coalescing compatible sequences where neither is subsumed in the
         * other but where the ancestral material is still entangled.
         */
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
            fprintf(print_progress, "Searching possible predecessors:\n");
        }
        // g_step_cost = 0;
        main_path_run_data.current_step_cost = 0;
        ac = COAL;
        preds = 0;
        nbdsize = 0;

        coalesce_compatibleandentangled_map(g, main_path_run_data, action);
        preds = predecessors.size() - nbdsize;
        nbdsize = predecessors.size();
        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Coalescing entangled: ", preds);
        }

        if (run_settings.se_cost != -1)
        {
            // g_step_cost = run_settings.se_cost;
            main_path_run_data.current_step_cost = run_settings.se_cost;
            ac = SE;

            seqerror_flips(g, main_path_run_data, action, run_settings);
            preds = predecessors.size() - nbdsize;
            nbdsize = predecessors.size();
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Sequencing errors: ", preds);
            }
        }

        if (run_settings.rm_cost != -1)
        {
            // g_step_cost = run_settings.rm_cost;
            main_path_run_data.current_step_cost = run_settings.rm_cost;
            ac = RM;

            recmut_flips(g, main_path_run_data, action, run_settings);

            preds = predecessors.size() - nbdsize;
            nbdsize = predecessors.size();
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Recurrent mutations: ", preds);
            }
        }

        /* Try all sensible events with one split */
        if (run_settings.r_cost != -1)
        {
            // g_step_cost = run_settings.r_cost;
            main_path_run_data.current_step_cost = run_settings.r_cost;
            ac = RECOMB1;

            maximal_prefix_coalesces_map(g, start, end, main_path_run_data, action);
            preds = predecessors.size() - nbdsize;
            nbdsize = predecessors.size();
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Prefix recombinations: ", preds);
            }

            maximal_postfix_coalesces_map(g, start, end, main_path_run_data, action);
            preds = predecessors.size() - nbdsize;
            nbdsize = predecessors.size();
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Postfix recombinations: ", preds);
            }
        }

        /* Try all sensible events with two splits */
        if (run_settings.rr_cost != -1)
        {
            // g_step_cost = run_settings.rr_cost;
            main_path_run_data.current_step_cost = run_settings.rr_cost;
            ac = RECOMB2;

            maximal_infix_coalesces_map(g, start, end, main_path_run_data, action);
            preds = predecessors.size() - nbdsize;
            nbdsize = predecessors.size();
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (infix): ", preds);
            }

            maximal_overlap_coalesces_map(g, start, end, main_path_run_data, action);
            preds = predecessors.size() - nbdsize;
            nbdsize = predecessors.size();
            if (g_howverbose > 0)
            {
                fprintf(print_progress, "%-40s %3d\n", "Two recombinations (overlap): ", preds);
            }
        }

        if (g_howverbose > 0)
        {
            fprintf(print_progress, "%-40s %3d\n", "Finished constructing predecessors.", predecessors.size());
            fprintf(print_progress, "-------------------------------------------------------------------------------------\n");
        }

        /* Finalise choice and prepare for next iteration */
        free_genes(g);

        /* Still looking for path to MRCA */
        if (!ontheflyselection && global)
        {
            /* So far we have only enumerated putative predecessors -
             * score these and choose one.
             */

            // Set the tracking lists to NULL for the score computation, and destroy the old g_sequence_labels/sites
            reset();

            nbdsize = predecessors.size(); // number of predecessors we score
            if (nbdsize == 0)
            {
                fprintf(stderr, "No neighbours left to search but MRCA not reached.");
            }
            total_nbdsize = total_nbdsize + nbdsize;

            // Calculate all the scores and store in an array
            // Update sc_min and sc_max for renormalising the score later
            score_array = (double *)malloc(predecessors.size() * sizeof(double));
            if (!_choice_fixed)
            {
                sc_min = DBL_MAX, sc_max = 0;
                int i = 0;
                for (auto &f : predecessors)
                {
                    RunData scoring_data(true); // Empty RunData that won't track anything
                    // output_genes(f->g, stderr, "\nGenes of f:\n");
                    reset_beagle_builtins(f->g, scoring_data); // set f to be _greedy_currentstate

                    // Calculate all the scores and update the min and max

                    score_array[i] = scoring_function(f->g, f->step_cost, run_settings, scoring_data);

                    i++;
                }
            }

            // Now consider each predecessor one by one, score, and set as the new choice if the score is lower
            int i = 0;
            for (auto &f : predecessors)
            {
                RunData scoring_data(true);
                reset_beagle_builtins(f->g, scoring_data); // set _greedy_currentstate to be f->g

                // g_step_cost = f->recombinations;
                printscore = score_renormalise(f->g, score_array[i], f->step_cost, run_settings, scoring_data);
                if (print_progress != NULL && g_howverbose == 2)
                {
                    fprintf(print_progress, "Predecessor %d obtained with event cost %.1f:\n", i + 1, f->step_cost);
                    output_genes(f->g, print_progress, NULL);
                    print_int_vector(f->elements, "Sequences: ");
                    print_int_vector(f->sites, "Sites: ");
                    fprintf(print_progress, "Predecessor score: %.0f \n\n",
                            (printscore == -DBL_MAX ? -INFINITY : (printscore == DBL_MAX ? INFINITY : printscore)));
                    fflush(print_progress);
                }
                if (select(printscore))
                {
                    // compute score and check if better than that of greedy_choice
                    /* Set f to be new choice */
                    greedy_choice = std::move(f);
                    // output_genes(greedy_choice->g, stderr, "greedy_choice update:\n");
                }

                i++;
            }

            free(score_array);
            predecessors.clear();
        }

        // greedy_choice is a unique pointer and everything it will be reset so need to copy out elements.
        // output_genes(greedy_choice->g, stderr, "greedy_choice:\n");
        g = greedy_choice->g;
        main_path_run_data.sequence_labels = greedy_choice->elements;
        main_path_run_data.site_labels = greedy_choice->sites;

        main_path_run_data.current_step_cost = greedy_choice->step_cost;

        switch (greedy_choice->action)
        {
        case COAL:
            break;
        case SE:
            seflips += greedy_choice->step_cost / run_settings.se_cost;
            break;
        case RM:
            rmflips += greedy_choice->step_cost / run_settings.rm_cost;
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

        if (g_use_eventlist && main_path_run_data.eventlist.in_use)
        {
            main_path_run_data.eventlist.append(greedy_choice->events);
        }

        r += greedy_choice->step_cost;

        /* Clean up */
        HistoryFragment *ptr = greedy_choice.release(); // Needed so that it doesn't delete all it's members (which other variables are now pointing to)
        free(ptr);
        free(start);
        free(end);
        if (_choice_fixed)
        {
            free_genes(g);
        }

        // Can abandon the run if the number of recombinations already exceeds rec_max
        if (recombs > run_settings.rec_max)
        {
            bad_soln = 1;
            break;
        }

        // Can also abandon the run if the number of SE+RM when we have r recombinations is greater than what we've
        // seen in earlier solutions.
        if (run_settings.rec_max != INT_MAX && !g_lookup.empty())
        {
            if (seflips + rmflips > g_lookup[recombs])
            {
                bad_soln = 1;
                break;
            }
        }
    }

    // If we exited the loop because of a sub-optimal solution, record this
    if (bad_soln)
    {
        if (run_settings.run_reference > 0)
        {
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", run_settings.run_reference, run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, total_nbdsize);
        }
        else
        {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f  NA  NA  NA %10d ", run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, total_nbdsize);
        }
    }
    else
    {
        // Otherwise, record the result
        if (print_progress != NULL && g_howverbose > 0)
        {
            fprintf(print_progress, "\nTotal number of states considered: %d\n", total_nbdsize);
            fprintf(print_progress, "Total event cost: %.1f\n", r);
            if (run_settings.run_reference > 0)
            {
                fprintf(print_progress, "%10s %13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Ref", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost",
                        "SE", "RM", "R", "N_states", "Time");
            }
            else
            {
                fprintf(print_progress, "%13s %6s %8s %8s %8s %8s %3s %3s %3s %10s %15s\n", "Seed", "Temp", "SE_cost", "RM_cost", "R_cost", "RR_cost", "SE", "RM", "R", "N_states", "Time");
            }
        }
        if (run_settings.run_reference > 0)
        {
            fprintf(print_progress, "%10d %13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", run_settings.run_reference, run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, seflips, rmflips, recombs, total_nbdsize);
        }
        else
        {
            fprintf(print_progress, "%13.0f %6.1f %8.2f %8.2f %8.2f %8.2f %3d %3d %3d %10d ", run_settings.run_seed, run_settings.temp, run_settings.se_cost, run_settings.rm_cost, run_settings.r_cost, run_settings.rr_cost, seflips, rmflips, recombs, total_nbdsize);
        }
        if (!g_lookup.empty())
        {
            if (seflips + rmflips < g_lookup[recombs])
            {
                // If found a better bound r < rec_max for Rmin, update.
                if (seflips + rmflips == 0)
                {
                    run_settings.rec_max = recombs;
                }
                update_lookup(g_lookup, recombs, seflips + rmflips);
            }
        }
    }

    return r;
}
