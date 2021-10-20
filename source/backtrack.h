/*******************************************************************

    backtrack.h

    Description of function and data types for backtracking a minimum
    recombination history

********************************************************************/

#ifndef BACKTRACK_H
#define BACKTRACK_H

#include "arg.h"
#include "gene.h"

typedef enum
{
    SUBSTITUTION,
    COALESCENCE,
    RECOMBINATION,
    REMOVE,
    COLLAPSE,
    SWAP,
    LOOKUP,
    SEFLIP,
    RMFLIP
} EventType;
typedef enum
{
    COAL,
    SE,
    RM,
    RECOMB1,
    RECOMB2
} Action;

typedef struct _Event
{
    EventType type;
    union
    {
        struct
        {
            int seq;
            int site;
        } s;
        struct
        {
            int s1;
            int s2;
        } c;
        struct
        {
            int seq;
            int pos;
        } r;
        struct
        {
            int s1;
            int s2;
        } swap;
        struct
        {
            int seq;
            int site;
        } flip;
        int collapse;
        int remove;
        int lookup;
    } event;
} Event;

typedef struct _HistoryFragment
{
    Genes *g;         /* End configuration */
    LList *event;     /* List of events leading from start
                       * configuration to end configuration.
                       */
    double step_cost; /* Number of recombination events */
    EList *elements;
    EList *sites;
    Action action;
} HistoryFragment;

typedef struct _PartialHistory
{
    Genes *g;               // current genes
    LList *event_list;      // List of Event's leading from start config to g
    EList *sequence_labels; // List of ints which label the sequences. starting ones are 0,1,... (become -1 as coalesce)
    EList *site_labels;     // List of ints which label the sites. starting ones are 0,1,... (become -1 as merged etc)
    int weight;             // used for doing dfs to keep track of how many leaves will come from here.
    int recombinations;
    int recurrent_mutations;
    int num_of_sequences;
} PartialHistory;

#ifdef DEBUG
extern HashTable *ancestral_state_trace;
#endif
ARG *eventlist2history(AnnotatedGenes *a, FILE *output);

#endif
