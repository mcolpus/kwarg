/*******************************************************************

    backtrack.h

    Description of function and data types for backtracking a minimum
    recombination history

********************************************************************/

#ifndef BACKTRACK_H
#define BACKTRACK_H

#include "arg.h"
#include "gene.h"

#include <vector>

typedef enum {SUBSTITUTION, COALESCENCE, RECOMBINATION, REMOVE,
    COLLAPSE, SWAP, LOOKUP, SEFLIP, RMFLIP} EventType;
typedef enum {COAL, SE, RM, RECOMB1, RECOMB2} Action;
    
    typedef struct _Event {
        EventType type;
        union {
            struct {
                int seq;
                int site;
            } s;
            struct {
                int s1;
                int s2;
            } c;
            struct {
                int seq;
                int pos;
            } r;
            struct {
                int s1;
                int s2;
            } swap;
            struct {
                int seq;
                int site;
            } flip;
            int collapse;
            int remove;
            int lookup;
        } event;
    } Event;
    

struct HistoryFragment {
  Genes *g;           /* End configuration */
  LList *event;       /* List of events leading from start
		       * configuration to end configuration.
		       */
  double recombinations; /* Number of recombination events */
  std::vector<int> elements;
  std::vector<int> sites;
  Action action;

  HistoryFragment() = default;

  HistoryFragment(HistoryFragment const &) = delete;
  HistoryFragment& operator=(HistoryFragment const &) = delete;

  ~HistoryFragment() {
    if(g != NULL)
      free_genes(g);
    if (event != NULL) {
      while (Length(event) != 0)
        free(Pop(event));
      DestroyLList(event);
    }
  }
};

#ifdef DEBUG
extern HashTable *ancestral_state_trace;
#endif
ARG *eventlist2history(AnnotatedGenes *a, FILE *output);

#endif
