/*******************************************************************

    backtrack.h

    Description of function and data types for backtracking a minimum
    recombination history

********************************************************************/

#ifndef BACKTRACK_H
#define BACKTRACK_H

#include <vector>
#include <list>

#include "arg.h"
#include "gene.h"

// Moved structs to gene to avoid cyclic references

#ifdef DEBUG
extern HashTable *ancestral_state_trace;
#endif
void output_g_eventlist_new_as_text(FILE *output);
ARG *eventlist2history(AnnotatedGenes *a, FILE *output);

#endif
