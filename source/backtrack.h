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
void output_eventlist_as_text(FILE *output, const EventList &eventlist);
ARG *eventlist2history(const AnnotatedGenes *a, FILE *output, EventList eventlist);

#endif
