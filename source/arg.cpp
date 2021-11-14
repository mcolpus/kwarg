/*******************************************************************
 *
 *    arg.c :Implementation of functions for building and outputting an
 *    ancestral recombination graph
 *
 ********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>

#include <list>
#include <vector>
#include <string>

#include "bitfunctions.h"
#include "llist.h"
#include "arg.h"
#include "common.h"


/* Add a type node to arg, with label and s as label and sequence. If
 * type is ARGRECOMBINATION, an extra argument specifying the first
 * site of the suffix contribution to the recombinant should be
 * provided. Return value is index of new node. The strings supplied
 * for label and sequence will be deallocated by arg_destroy.
 */
int arg_addnode(ARG &arg, ARGNodeType type, std::string label, std::string seq, int pos)
{
    va_list args;

    ARGNode new_node;

    new_node.type = type;
    new_node.label = label;
    new_node.sequence = seq;

    switch (type)
    {
    case ARGSAMPLE:
    case ARGCOALESCENCE:
        /* These types of nodes have one predecessor */
        new_node.predecessor.one.mutations = {};
        new_node.predecessor.one.target = -1;
        break;
    case ARGRECOMBINATION:
        /* These types of nodes have two predecessors */
        new_node.predecessor.two.prefix.mutations = {};
        new_node.predecessor.two.prefix.target = -1;
        new_node.predecessor.two.suffix.target = -1;
        new_node.predecessor.two.position = pos;
        new_node.predecessor.two.suffix.mutations = {};
        break;
    case ARGANCESTOR:
        /* These types of nodes have no predecessors */
        break;
    }

    arg.nodes.push_back(new_node);

    return arg.size();
}

int arg_addnode(ARG &arg, ARGNodeType type, std::string label, std::string seq)
{
    return arg_addnode(arg, type, label, seq, 0);
}

/* Return node with index i in arg */
ARGNode *arg_getnode(ARG &arg, int i)
{
    return &(arg.nodes[i]);
}

/* Finalise arg, changing nodes with non-terminated outgoing edges to
 * ANCESTOR type.
 */
void arg_finalise(ARG &arg)
{
    int i;

    for (i = 0; i < arg.size(); i++)
    {
        switch (arg.nodes[i].type)
        {
        case ARGSAMPLE:
        case ARGCOALESCENCE:
            if (arg.nodes[i].predecessor.one.target == -1)
            {
                /* Sanity check */
                if (arg.nodes[i].predecessor.one.mutations.size() > 0)
                {
                    fprintf(stderr, "Something is wrong in ARG finalisation - please email error report\n");
                    exit(1);
                }
                arg.nodes[i].type = ARGANCESTOR;
            }
            break;
        case ARGRECOMBINATION:
            if (arg.nodes[i].predecessor.two.prefix.target == -1)
                arg.nodes[i].predecessor.two.prefix.target = arg_addnode(arg, ARGANCESTOR, "NULL", "NULL");
            if (arg.nodes[i].predecessor.two.suffix.target == -1)
                arg.nodes[i].predecessor.two.suffix.target = arg_addnode(arg, ARGANCESTOR, "NULL", "NULL");
            break;
        case ARGANCESTOR:
            break;
        }
    }
}

/* Output labels to fp, separated by semicolons */
static void output_edgelabels(std::list<int> &labels, AnnotatedGenes *a, bool maiden,
                              FILE *fp)
{
    for (const int p : labels)
    {
        if (!maiden)
            fprintf(fp, ";");
        if (p == -INT_MAX)
        {
            if (a->positions != NULL)
                fprintf(fp, "*%s", (char *)GetByIndex(a->positions, 0));
            else
                fprintf(fp, "*%d", 1);
        }
        else
        {
            if (p < 0)
            {
                if (a->positions != NULL)
                    fprintf(fp, "*%s", (char *)GetByIndex(a->positions, -p));
                else
                    fprintf(fp, "*%d", -p + 1);
            }
            else
            {
                if (a->positions != NULL)
                    fprintf(fp, "%s", (char *)GetByIndex(a->positions, p));
                else
                    fprintf(fp, "%d", p + 1);
            }
        }
        maiden = false;
    }
}

static void update_descendant(int (*a)[2], int d)
{
    if ((*a)[0] == -1)
        (*a)[0] = d;
    else
        (*a)[1] = d;
}

static void transfer_mutations(std::list<int> &from, std::list<int> &to, int startsite,
                               int endsite)
{
    for (const int site : from)
    {
        if ((site >= startsite) && (site < endsite))
            to.push_back(site);
    }
}

static void transfer_mutations(LList *from, LList *to, int startsite,
                               int endsite)
{
    LListCounter *lcounter = MakeCounter(from, FIRST);
    int i, site;

    for (i = 0; i < Length(from); i++)
    {
        site = (intptr_t)Next(lcounter);
        if ((site >= startsite) && (site < endsite))
            Enqueue(to, (void *)site);
    }
    free(lcounter);
}

static void visit(ARG &arg, ARG &tree, int current, int (*descendant)[2],
                  int *tree_node, int startsite, int endsite)
{
    switch (arg.nodes[current].type)
    {
    case ARGSAMPLE:
        update_descendant(descendant + arg.nodes[current].predecessor.one.target,
                          current);
        transfer_mutations(arg.nodes[current].predecessor.one.mutations,
                           tree.nodes[tree_node[current]].predecessor.one.mutations,
                           startsite, endsite);
        visit(arg, tree, arg.nodes[current].predecessor.one.target, descendant,
              tree_node, startsite, endsite);
        break;
    case ARGCOALESCENCE:
    case ARGANCESTOR:
        if (descendant[current][1] == -1)
            /* We still need to visit the second descendant of this node */
            return;
        if (tree_node[descendant[current][0]] != -1)
        {
            if (tree_node[descendant[current][1]] != -1)
            {
                /* Both descendants are part of tree */
                tree.nodes.push_back(arg.nodes[current]);
                // memcpy(tree.nodes + tree.size(), arg.nodes + current, sizeof(ARGNode));
                tree.nodes[tree_node[descendant[current][0]]].predecessor.one.target = tree.size();
                tree.nodes[tree_node[descendant[current][1]]].predecessor.one.target = tree.size();
                if (arg.nodes[current].type != ARGANCESTOR)
                {
                    tree.nodes[tree.size()].predecessor.one.mutations = {};
                    transfer_mutations(arg.nodes[current].predecessor.one.mutations,
                                       tree.nodes[tree.size()].predecessor.one.mutations,
                                       startsite, endsite);
                    update_descendant(descendant + arg.nodes[current].predecessor.one.target, current);
                }
                tree_node[current] = tree.size();
            }
            else
            {
                /* Only first descendant is part of tree */
                if (arg.nodes[current].type == ARGANCESTOR)
                {
                    /* First descendant is common ancestor of all samples */
                    tree.nodes[tree_node[descendant[current][0]]].type = ARGANCESTOR;
                }
                else
                {
                    /* Descendant of ancestor to current node in tree is first
                     * descendant to current node.
                     */
                    update_descendant(descendant + arg.nodes[current].predecessor.one.target,
                                      descendant[current][0]);
                    /* Copy relevant mutations from edge to predecessor in ARG
                     * to edge from first descendant in marginal tree.
                     */
                    transfer_mutations(arg.nodes[current].predecessor.one.mutations,
                                       tree.nodes[tree_node[descendant[current][0]]].predecessor.one.mutations,
                                       startsite, endsite);
                }
            }
        }
        else
        {
            if (tree_node[descendant[current][1]] != -1)
            {
                /* Only second descendant is part of tree */
                if (arg.nodes[current].type == ARGANCESTOR)
                {
                    /* Second descendant is common ancestor of all samples */
                    tree.nodes[tree_node[descendant[current][1]]].type = ARGANCESTOR;
                }
                else
                {
                    /* Descendant of ancestor to current node in tree is second
                     * descendant to current node.
                     */
                    update_descendant(descendant + arg.nodes[current].predecessor.one.target,
                                      descendant[current][1]);
                    /* Copy relevant mutations from edge to predecessor in ARG
                     * to edge from first descendant in marginal tree.
                     */
                    transfer_mutations(arg.nodes[current].predecessor.one.mutations,
                                       tree.nodes[tree_node[descendant[current][1]]].predecessor.one.mutations,
                                       startsite, endsite);
                }
            }
            else
                /* None of the descendants are part of the tree */
                update_descendant(descendant + arg.nodes[current].predecessor.one.target, current);
        }
        if (arg.nodes[current].type != ARGANCESTOR)
            visit(arg, tree, arg.nodes[current].predecessor.one.target, descendant,
                  tree_node, startsite, endsite);
        break;
    case ARGRECOMBINATION:
        if (startsite < arg.nodes[current].predecessor.two.position)
        {
            /* Predecessor providing prefix is predecessor at this site */
            update_descendant(descendant + arg.nodes[current].predecessor.two.prefix.target,
                              descendant[current][0]);
            /* Copy relevant mutations from edge to predecessor in ARG
             * to edge from first descendant in marginal tree.
             */
            if (tree_node[descendant[current][0]] != -1)
            {
                transfer_mutations(arg.nodes[current].predecessor.two.prefix.mutations,
                                   tree.nodes[tree_node[descendant[current][0]]].predecessor.one.mutations,
                                   startsite, endsite);
            }

            visit(arg, tree, arg.nodes[current].predecessor.two.prefix.target,
                  descendant, tree_node, startsite, endsite);
            /* Visit non-predecessor */
            update_descendant(descendant + arg.nodes[current].predecessor.two.suffix.target,
                              current);
            visit(arg, tree, arg.nodes[current].predecessor.two.suffix.target,
                  descendant, tree_node, startsite, endsite);
        }
        else
        {
            /* Predecessor providing suffix is predecessor at this site */
            update_descendant(descendant + arg.nodes[current].predecessor.two.suffix.target,
                              descendant[current][0]);
            /* Copy relevant mutations from edge to predecessor in ARG
             * to edge from first descendant in marginal tree.
             */
            if (tree_node[descendant[current][0]] != -1)
                transfer_mutations(arg.nodes[current].predecessor.two.suffix.mutations,
                                   tree.nodes[tree_node[descendant[current][0]]].predecessor.one.mutations,
                                   startsite, endsite);
            visit(arg, tree, arg.nodes[current].predecessor.two.suffix.target,
                  descendant, tree_node, startsite, endsite);
            /* Visit non-predecessor */
            update_descendant(descendant + arg.nodes[current].predecessor.two.prefix.target,
                              current);
            visit(arg, tree, arg.nodes[current].predecessor.two.prefix.target,
                  descendant, tree_node, startsite, endsite);
        }
        break;
    }
}

static void check_positions(int site, int *quote, LList *labels)
{
    char *label = (char *)GetByIndex(labels, site);

    /* If we have already seen a single quote there is no point in
     * looking for more.
     */
    if (!*quote)
        for (int i = 0; i < strlen(label); i++)
            if (label[i] == '\'')
                *quote = 1;
}

static void newick_visit(ARG &tree, int node, AnnotatedGenes *a,
                         ARGLabels nodelabels, int annotate_edges,
                         int generate_id, std::string parts[], std::string subtree)
{
    int i, j, single_quotes;
    std::string label = tree.nodes[node].label;
    LListCounter *lcounter;

    if ((tree.nodes[node].type == ARGSAMPLE) || (parts[node] != "NULL"))
    {
        /* Last visit to this node */
        if (tree.nodes[node].type != ARGSAMPLE)
        {
            /* Finish construction of bifucating node */
            parts[node] += ',';
            if (annotate_edges)
                /* Add space between nodes to improve readability */
                parts[node] += ' ';
            parts[node] += subtree;
            parts[node] += ')';
        }

        /* Add label to node */
        if (((!label.empty()) || ((tree.nodes[node].type == ARGSAMPLE) && (generate_id))) && (nodelabels != ARGNONE) && (nodelabels != ARGSEQUENCE) && ((nodelabels != ARGSEQUENCEFIRST) || (tree.nodes[node].sequence == "NULL")))
        {
            /* The label of the node is used for actual labelling */
            /* If label contains single quotes we need to escape them */
            single_quotes = 0;
            if (!label.empty())
            {
                /* Genuine label */
                for (i = 0; label[i] != '\0'; i++)
                    if (label[i] == '\'')
                        single_quotes++;
            }
            else
                /* Label is just sample index */
                label = i2a(node + 1);

            if (single_quotes > 0)
            {
                /* Label did contain single quotes */
                j = single_quotes;
                label = (char *)xmalloc((i + j + 1) * sizeof(char));
                for (; i + j >= 0; i--)
                {
                    label[i + j] = label[i];
                    if (label[i] == '\'')
                        label[i + --j] = '\'';
                }
            }

            /* Add label annotation to Newick representation */
            if (((nodelabels == ARGBOTH) && (!tree.nodes[node].sequence.empty())) || (single_quotes > 0))
                /* Node label needs to be quoted */
                parts[node] += '\'';
            parts[node] += tree.nodes[node].label;
            if ((nodelabels == ARGBOTH) && (!tree.nodes[node].sequence.empty()))
            {
                /* Both label and sequence should be used as annotation */
                parts[node] += ';' + tree.nodes[node].sequence + '\'';
            }
            else if (single_quotes > 0)
                /* End label quote */
                parts[node] += '\'';

        }
        else if ((tree.nodes[node].sequence != "NULL") &&
                 (nodelabels != ARGNONE) && (nodelabels != ARGLABEL))
        {
            /* The sequence of the node is used for labelling */
            parts[node] += tree.nodes[node].sequence;
        }

        if (tree.nodes[node].type != ARGANCESTOR)
        {
            if (annotate_edges)
            {
                /* Add mutations on edge from node to predecessor */
                if (tree.nodes[node].predecessor.one.mutations.size() > 0)
                {
                    parts[node] += ':';

                    /* Check whether list of mutations needs to be quoted */
                    single_quotes = 0;
                    if (a->positions != NULL)
                    {
                        for (const int mut : tree.nodes[node].predecessor.one.mutations)
                        {
                            check_positions(mut, &single_quotes, a->positions);
                        }
                    }

                    if (single_quotes)
                        /* It does */
                        parts[node] += '\'';
                    /* Add each mutation to 'length' of this node's ancestral edge */

                    bool first_loop = true;
                    for (const int j : tree.nodes[node].predecessor.one.mutations)
                    {
                        if (first_loop)
                        {
                            parts[node] += ';';
                            first_loop = false;
                        }

                        /* Get next mutation and add it to length label */
                        if (a->positions != NULL)
                        {
                            if (single_quotes)
                            {
                                /* Quote the single quotes in each site label */
                                label = (char *)GetByIndex(a->positions, j);
                                for (int i = 0; i < label.size(); i++)
                                {
                                    parts[node] += label[i];
                                    if (label[i] == '\'')
                                        parts[node] += '\'';
                                }
                            }
                            else
                                parts[node] += (char *)GetByIndex(a->positions, j);
                        }
                        else
                        {
                            label = i2a(j + 1);
                            parts[node] += label;
                        }
                    }
                }
            }

            /* Recurse on parent node */
            if (tree.nodes[node].predecessor.one.target < tree.size())
            {
                newick_visit(tree, tree.nodes[node].predecessor.one.target, a, nodelabels, annotate_edges, generate_id, parts, parts[node]);
            }
        }
    }
    else
    {
        /* First visit to bifucating node - initialise construction */
        parts[node] = "(";
        parts[node] += subtree;
    }
}

static void newick_output(ARG &tree, AnnotatedGenes *a, FILE *fp,
                          ARGLabels nodelabels, int annotate_edges,
                          int generate_id, char *label)
{
    int i;
    std::string parts[tree.size()];
    //MyString **parts = new MyString *[tree.size()];
    std::string s;

    /* Initialise parts */
    for (i = 0; i < tree.size(); i++)
        if (tree.nodes[i].type == ARGSAMPLE)
            parts[i] = "";
        else
            parts[i] = "NULL";

    /* Generate newick representation */
    for (i = 0; i < tree.size(); i++)
        if (tree.nodes[i].type == ARGSAMPLE)
            newick_visit(tree, i, a, nodelabels, annotate_edges, generate_id, parts,
                         NULL);

    /* Output Newick representation */
    /* Ancestor should always be last node in ARG structure */
    s = parts[tree.size() - 1];
    fputs(s.c_str(), fp);
    if (label != NULL)
        fputs(label, fp);
    putc(';', fp);
    putc('\n', fp);

    /* Clean up */
    /* All MyStrings in parts should have been destroyed by
     * concatenations to other MyStrings, except for the one at the
     * ancestor.
     */
}

static void marginal_trees(ARG &arg, AnnotatedGenes *a, FILE *fp,
                           ARGFormat format, ARGLabels nodelabels,
                           int annotate_edges, int generate_id, int intervals,
                           char *label)
{
    int i, j, next, *starts, *tree_node, (*descendant)[2];
    char *t;
    MyString *s;

    if (intervals)
    {
        /* Start by determining the recombination free intervals */
        starts = (int *)xcalloc(a->g->length, sizeof(int));
        starts[0] = 1;
        for (i = 0; i < arg.size(); i++)
            if (arg.nodes[i].type == ARGRECOMBINATION)
                starts[arg.nodes[i].predecessor.two.position] = 1;
        /* Convert Booleans to index of next interval start */
        next = a->g->length;
        for (i = a->g->length - 1; i >= 0; i--)
            if (starts[i])
            {
                starts[i] = next;
                next = i;
            }
    }
    else
    {
        /* Each position should be treated as its own interval */
        starts = (int *)xmalloc(a->g->length * sizeof(int));
        for (i = 0; i < a->g->length; i++)
            starts[i] = i + 1;
    }

    /* Generate and output evolutionary tree(s) for each recombination
     * free interval.
     */
    ARG tree;
    i = 0;
    tree_node = (int *)xmalloc(arg.size() * sizeof(int));
    descendant = (int(*)[2])xmalloc(arg.size() * sizeof(int[2]));
    /* We need to copy all sample nodes to the marginal tree first to
     * have generate_id have the desired effect.
     */
    for (j = 0; j < arg.size(); j++)
        if (arg.nodes[j].type == ARGSAMPLE)
        {
            tree.nodes.push_back(arg.nodes[j]);
            // memcpy(tree.nodes + tree.size(), arg.nodes + j, sizeof(ARGNode));
            tree.nodes[tree.size()].predecessor.one.mutations = {};
            tree_node[j] = tree.size();
        }

    /* Continue until we reach the end of the sequences */
    while (i < a->g->length)
    {
        /* Wipe all non-sample nodes from tree */
        // tree.n = a->g->n;
        /* Initialise arrays of descendants of a node and for mapping from
         * ARG nodes to marginal tree nodes.
         */
        for (j = 0; j < arg.size(); j++)
            if (arg.nodes[j].type != ARGSAMPLE)
            {
                descendant[j][0] = descendant[j][1] = -1;
                tree_node[j] = -1;
            }

        /* Start traversal back to root from each sample node */
        for (j = 0; j < arg.size(); j++)
            if (arg.nodes[j].type == ARGSAMPLE)
                visit(arg, tree, j, descendant, tree_node, i, starts[i]);

        /* Copy sequence to tree nodes */
        for (j = 0; j < arg.size(); j++)
            if ((tree_node[j] != -1) && (arg.nodes[j].sequence != "NULL"))
            {
                /* Copy relevant part of the sequence from ARG node to tree node */
                tree.nodes[tree_node[j]].sequence = (char *)xmalloc((starts[i] - i + 1) * sizeof(char));
                tree.nodes[tree_node[j]].sequence = arg.nodes[j].sequence;
            }
        /* Output current tree */
        /* GDL: graph can vaere et element i en anden graph; title: text angiver
         *      titel paa graf
         * GML: */
        /* Generate label specifying site or region tree is valid for */
        if (a->positions != NULL)
        {
            s = mystring_str2mystr((char *)GetByIndex(a->positions, i));
            if (i < starts[i] - 1)
            {
                /* Tree valid for interval rather than single site */
                mystring_addend(s, '-');
                mystring_append(s, (char *)GetByIndex(a->positions, starts[i] - 1));
            }
        }
        else
        {
            t = i2a(i + 1);
            s = mystring_str2mystr(t);
            free(t);
            if (i < starts[i] - 1)
            {
                /* Tree valid for interval rather than single site */
                mystring_addend(s, '-');
                t = i2a(starts[i]);
                mystring_append(s, t);
                free(t);
            }
        }
        t = mystring_mystr2str(s);
        mystring_free(s);
        if (format == TREENEWICK)
            newick_output(tree, a, fp, nodelabels, annotate_edges, generate_id, t);
        else if (format == TREEDOT)
            arg_output(tree, a, fp, ARGDOT, nodelabels, annotate_edges, generate_id);
        else if (format == TREEGDL)
            arg_output(tree, a, fp, ARGGDL, nodelabels, annotate_edges, generate_id);
        else if (format == TREEGML)
            arg_output(tree, a, fp, ARGGML, nodelabels, annotate_edges, generate_id);
        /* Clean up */
        free(t);
        /* Clear everything but non-site specific information of sample
         * nodes from marginal tree data structure.
         */
        for (j = 0; j < tree.size(); j++)
        {
            if (tree.nodes[j].type != ARGANCESTOR)
            {
                tree.nodes[j].predecessor.one.mutations.clear();
            }
        }
        i = starts[i];
    }

    /* Clean up */
    for (j = 0; j < tree.size(); j++)
        if (tree.nodes[j].type == ARGSAMPLE)
            tree.nodes[j].predecessor.one.mutations.clear();
    free(starts);
    free(tree_node);
    free(descendant);
}

/* Output arg to file fp (stdout if fp is NULL). Sequence and site
 * information is taken from a, graph description language (DOT, GDL
 * or GML) is specified by format, node labelling (see definition of
 * ARGLabels in arg.h) is specified by nodelabels; if annotate_edges
 * is true edges are labelled with mutations, and if generate_id is
 * true sample nodes with no label in a are assigned their index in a
 * as label. If format is one of the formats prefixed with TREE an
 * extra argument should be provided, specifying whether one tree
 * for every site (if argument is false) or one tree for every
 * recombination free interval (if argument is true) should be output.
 */
void arg_output(ARG &arg, AnnotatedGenes *a, FILE *fp, ARGFormat format,
                ARGLabels nodelabels, int annotate_edges, int generate_id, ...)
{
    int i, maiden;
    va_list args;
    int intervals;

    /* Check whether marginal trees rather than the arg should be output */
    switch (format)
    {
    case TREEDOT:
    case TREEGDL:
    case TREEGML:
    case TREENEWICK:
        /* Optional argument should specify whether a tree should be
         * output for every site, or only for every recombination free
         * interval.
         */
        va_start(args, generate_id);
        intervals = (int)va_arg(args, int);
        va_end(args);
        marginal_trees(arg, a, fp, format, nodelabels, annotate_edges, generate_id,
                       intervals, NULL);
        return;
    }

    if (fp == NULL)
        fp = stdout;

    /* Output prelude */
    switch (format)
    {
    case ARGDOT:
        fprintf(fp, "digraph ARG {\n  { rank = same;");
        for (i = 0; i < a->g->n; i++)
            fprintf(fp, " %d;", i);
        fprintf(fp, " }\n");
        break;
    case ARGGDL:
        fprintf(fp, "graph: {\n");
        break;
    case ARGGML:
        fprintf(fp, "graph [\n  directed 1\n");
        break;
    }
    /* Output nodes and their edges */
    for (i = 0; i < arg.size(); i++)
    {
        /* Output node i */
        maiden = 1;
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, "  %d [label=\"", i);
            break;
        case ARGGDL:
            fprintf(fp, "  node: { title: \"%d\" label: \"", i);
            break;
        case ARGGML:
            fprintf(fp, "  node [\n    id %d\n    label \"", i);
            break;
        }
        /* Generate node label */
        if (((arg.nodes[i].label != "NULL") || ((arg.nodes[i].type == ARGSAMPLE) && generate_id)) && ((nodelabels == ARGLABEL) || (nodelabels == ARGBOTH) || (nodelabels == ARGLABELFIRST) || ((nodelabels == ARGSEQUENCEFIRST) && (arg.nodes[i].sequence == "NULL"))))
        {
            if (arg.nodes[i].label != "NULL")
                fprintf(fp, "%s", arg.nodes[i].label);
            else
                fprintf(fp, "%d", i + 1);
            maiden = 0;
        }
        if ((arg.nodes[i].sequence != "NULL") && ((nodelabels == ARGSEQUENCE) || (nodelabels == ARGBOTH) || (nodelabels == ARGSEQUENCEFIRST) || ((nodelabels == ARGLABELFIRST) && maiden)))
        {
            if (!maiden)
                fprintf(fp, "; %s", arg.nodes[i].sequence);
            else
                fprintf(fp, "%s", arg.nodes[i].sequence);
            maiden = 0;
        }
        if (maiden)
        {
            /* No node label printed */
            if ((arg.nodes[i].type == ARGRECOMBINATION) && (arg.nodes[i].label != "NULL"))
                /* If a recombination node is still unlabeled, label it (label
                 * should be recombination point) if possible.
                 */
                fprintf(fp, "%s\"", arg.nodes[i].label);
            else
            {
                switch (format)
                {
                case ARGDOT:
                    fprintf(fp, "\",shape=point");
                    break;
                case ARGGDL:
                    fprintf(fp, "\" scaling: 0.0");
                    break;
                case ARGGML:
                    fprintf(fp, " \"");
                    break;
                }
            }
        }
        else
            fprintf(fp, "\"");
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, ",color=");
            break;
        case ARGGDL:
            fprintf(fp, " shape: circle bordercolor: ");
            break;
        case ARGGML:
            fprintf(fp, "\n    graphics [\n      outline \"");
            break;
        }
        switch (arg.nodes[i].type)
        {
        case ARGSAMPLE:
            fprintf(fp, "red");
            break;
        case ARGCOALESCENCE:
            fprintf(fp, "green");
            break;
        case ARGRECOMBINATION:
            fprintf(fp, "blue");
            break;
        case ARGANCESTOR:
            fprintf(fp, "yellow");
            break;
        }
        /* End node description */
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, "];\n");
            break;
        case ARGGDL:
            if (i < a->g->n)
                fprintf(fp, " vertical_order: maxlevel");
            fprintf(fp, " }\n");
            break;
        case ARGGML:
            fprintf(fp,
                    "\"\n    ]\n    vgj [\n      type \"Oval\"\n      labelPosition \"in\"\n    ]\n  ]\n");
            break;
        }
        /* Output edges out of this node */
        switch (arg.nodes[i].type)
        {
        case ARGSAMPLE:
        case ARGCOALESCENCE:
            /* One outgoing edge */
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d", arg.nodes[i].predecessor.one.target, i);
                break;
            case ARGGDL:
                fprintf(fp, "  edge: { sourcename: \"%d\" targetname: \"%d\"",
                        arg.nodes[i].predecessor.one.target, i);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d",
                        arg.nodes[i].predecessor.one.target, i);
                break;
            }
            /* Output edge label */
            if (annotate_edges && arg.nodes[i].predecessor.one.mutations.size() > 0)
            {
                switch (format)
                {
                case ARGDOT:
                    fprintf(fp, " [label=\"");
                    break;
                case ARGGDL:
                    fprintf(fp, " label: \"");
                    break;
                case ARGGML:
                    fprintf(fp, "\n    label \"");
                    break;
                }
                output_edgelabels(arg.nodes[i].predecessor.one.mutations, a, true, fp);
                fprintf(fp, "\"");
                if (format == ARGDOT)
                    fprintf(fp, "]");
            }
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, ";\n");
                break;
            case ARGGDL:
                fprintf(fp, " }\n");
                break;
            case ARGGML:
                fprintf(fp, "\n  ]\n");
                break;
            }
            break;
        case ARGRECOMBINATION:
            /* Two outgoing edges */
            /* Prefix edge */
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d [label=\"P",
                        arg.nodes[i].predecessor.two.prefix.target, i);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"P",
                        arg.nodes[i].predecessor.two.prefix.target, i);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"P",
                        arg.nodes[i].predecessor.two.prefix.target, i);
                break;
            }

            if (annotate_edges && arg.nodes[i].predecessor.two.prefix.mutations.size() > 0)
                output_edgelabels(arg.nodes[i].predecessor.two.prefix.mutations, a, false, fp);

            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "\"]\n");
                break;
            case ARGGDL:
                fprintf(fp, "\" }\n");
                break;
            case ARGGML:
                fprintf(fp, "\"\n  ]\n");
                break;
            }
            /* Suffix edge */
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d [label=\"S",
                        arg.nodes[i].predecessor.two.suffix.target, i);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"S",
                        arg.nodes[i].predecessor.two.suffix.target, i);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"S",
                        arg.nodes[i].predecessor.two.suffix.target, i);
                break;
            }

            if (annotate_edges && arg.nodes[i].predecessor.two.suffix.mutations.size() > 0)
                output_edgelabels(arg.nodes[i].predecessor.two.suffix.mutations, a, false, fp);

            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "\"]\n");
                break;
            case ARGGDL:
                fprintf(fp, "\" }\n");
                break;
            case ARGGML:
                fprintf(fp, "\"\n  ]\n");
                break;
            }
            break;
        case ARGANCESTOR:
            /* No outgoing edges */
            break;
        }
    }
    /* Output postlude */
    switch (format)
    {
    case ARGDOT:
    case ARGGDL:
        fprintf(fp, "}\n");
        break;
    case ARGGML:
        fprintf(fp, "]\n");
        break;
    }
}
