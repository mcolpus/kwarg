#ifndef GENE_H
#define GENE_H

#include <stdio.h>

#include <vector>
#include <functional>
#include <memory>
#include <list>

#include "llist.h"
#include "hashtable.h"

/* Datatypes */
typedef struct _Gene {
  unsigned long *type;      /* Type bit vector */
  unsigned long *ancestral; /* Ancestral site bit vector */
} Gene;

typedef struct _Genes {
  int n;      /* Number of sequences */
  int length; /* Length of sequences */
  Gene *data; /* Sequences */
} Genes;

typedef struct _AnnotatedGenes {
  Genes *g;         /* Data */
  LList *positions; /* Site labels */
  LList *sequences; /* Sequence labels */
} AnnotatedGenes;

typedef enum { GENE_ANY, GENE_BEAGLE, GENE_FASTA } Gene_Format;
typedef enum { GENE_BINARY, GENE_NUCLEIC, GENE_AMINO } Gene_SeqType;

typedef struct _PackedGenes {
  int n;               /* Number of sequences */
  int length;          /* Length of sequences */
  unsigned int *data; /* Sequences */
} PackedGenes;

typedef struct _Site {
  unsigned long *type;      /* Type bit vector */
  unsigned long *ancestral; /* Ancestral sequence bit vector */
} Site;

typedef struct _Sites {
  int n;      /* Number of sequences */
  int length; /* Length of sequences */
  Site *data; /* Sites */
} Sites;

typedef struct _Index {
  int index;
  int block;
} Index;

/* Global variable for specifying whether the common ancestral
 * sequence is known or not.
 */
extern int gene_knownancestor;
/* Prototypes */
void free_genes(Genes *g);

// Moved from backtrack to avoid cyclic references
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

typedef struct _EVENTLIST
{
    bool in_use; // set to false when not wanting to record some section of computation
    std::list<Event> events;

    _EVENTLIST()
    {
        in_use = true;
        events = {};
    }

    int size()
    {
        if(!in_use)
        {
            fprintf(stderr, "EventList is not in use!");
        }
        return events.size();
    }

    void push_back(const Event &e)
    {
        if(!in_use)
        {
            fprintf(stderr, "EventList is not in use!");
        }
        events.push_back(e);
    }

    void set_null()
    {
        events = {};
        in_use = false;
    }

    void reset()
    {
        events = {};
        in_use = true;
    }

    // splices elements from other list to front
    void prepend(struct _EVENTLIST &other)
    {
        events.splice(events.begin(), other.events);
    }

    // splices elements from other list to front
    void append(struct _EVENTLIST &other)
    {
        events.splice(events.end(), other.events);
    }

    void destroy()
    {
        in_use = false;
    }

} EventList;

struct HistoryFragment {
  Genes *g;           /* End configuration */
  EventList events;       /* List of events leading from start
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
    events.destroy();
  }
};

#include "backtrack.h"

void free_annotatedgenes(AnnotatedGenes *g);
void free_sites(Sites *s);
void free_gene(Gene *g);
void free_site(Site *s);
Genes *new_genes(int n, int length, char **s, Gene_SeqType t);
AnnotatedGenes *read_genes(char *fname, Gene_Format f, Gene_SeqType t);
void output_genes(Genes *g, FILE *fp, char *comment);
void output_labelled_genes(Genes *g, FILE *fp, LList *labels);
void output_genes_indexed(Genes *s, FILE *fp);
void output_annotatedgenes(AnnotatedGenes *a, FILE *fp, char *comment);
char get_genes_character(Genes *g, int seq, int site);
void set_genes_character(Genes *g, int seq, int site, char c);
void swap_genes(Genes *g, int a, int b);
void add_ancestral_sites(Sites *s);
Genes *copy_genes(Genes *g);
Sites *copy_sites(Sites *s);
Genes *copy_allbutone(Genes *g, int a);
Genes *copy_region(Genes *g, int a, int b);
void remove_gene(Genes *g, int a);
void remove_annotatedgene(AnnotatedGenes *g, int a);
void reallocate_genes(Genes *g);
Sites *genes2sites(Genes *g);
char **genes2string(Genes *g);
void output_sites(Sites *s, FILE *fp, char *comment);
int remove_siamesetwins(Genes *g);
int remove_uninformative(Genes *g);
int remove_nonsegregating(Genes *g);
int coalesce_subsumed(Genes *g);
int implode_genes(Genes *g);
int no_recombinations_required(Genes *g);
void force_safeevents(Genes *g);
int force_mutations(Genes *g);
int mutate(Genes *g, int pos, int mutant);
std::vector<Genes *> force_mutation(Genes *g, std::vector<Event> &events);
std::vector<Genes *> force_mutation(Genes *g);
int segregating_site(Genes *g, int i);
int compatible(Genes *g, int a, int b);
std::vector<int> incompatible_sites(Genes *g, int a, int b);
int find_safe_coalescence(Genes *g, int a);
int entangled(Genes *g, int a, int b);
void coalesce(Genes *g, int a, int b);
std::vector<Genes *> force_coalesce(Genes *g, std::vector<Event> &events);
void split(Genes *g, int a, int i);
std::vector<Genes *> force_split(Genes *g, int a);
std::vector<Genes *> force_split(Genes *g, int a, std::vector<Event> &events);
void split_coalesceprefix(Genes *g, int a, int index, int block, int b);
void split_coalescepostfix(Genes *g, int a, int index, int block, int b);
void splitafter_coalescepostfix(Genes *g, int a, int index, int block, int b);
int split_removeprefix(Genes *g, int a, int index, int block);
int split_removepostfix(Genes *g, int a, int index, int block);
int ancestral_material(Genes *g);
int individual_all_ancestral(Genes *g, int i);
int ancestral_material_overlap(Genes *g);
int minimum_compatiblechops(Genes *g, int a);
Index *maximumsubsumedprefixs(Genes *g);
Index *maximumsubsumedpostfixs(Genes *g);
Index *maximumsubsumedprefix(Genes *g, int s);
Index *maximumsubsumedpostfix(Genes *g, int s);
std::vector<std::unique_ptr<HistoryFragment>> maximal_prefix_coalesces(Genes *g, Index *a, Index *b);
void maximal_prefix_coalesces_map(Genes *g, Index *a, Index *b,
				    std::function<void (Genes *)> f);
std::vector<std::unique_ptr<HistoryFragment>> maximal_postfix_coalesces(Genes *g, Index *a, Index *b);
void maximal_postfix_coalesces_map(Genes *g, Index *a, Index *b,
				     std::function<void (Genes *)> f);
void seqerror_flips(Genes* g, std::function<void (Genes *)> f);
void recmut_flips(Genes* g, std::function<void (Genes *)> f);
std::vector<std::unique_ptr<HistoryFragment>> maximal_infix_coalesces(Genes *g, Index *a, Index *b);
void maximal_infix_coalesces_map(Genes *g, Index *a, Index *b,
				  std::function<void (Genes *)> f);
std::vector<std::unique_ptr<HistoryFragment>> maximal_overlap_coalesces(Genes *g, Index *a, Index *b);
void maximal_overlap_coalesces_map(Genes *g, Index *a, Index *b,
				     std::function<void (Genes *)>);
int compare_sites(Sites *s, int a, int b);
int compare_genes(Genes *g, Genes *h);
PackedGenes *pack_genes(Genes *g);
void free_packedgenes(PackedGenes *p);
Genes *unpack_genes(PackedGenes *p);
int compare_packedgenes(PackedGenes *g, PackedGenes *h);
HashTable *new_packedgeneshashtable(int bits);



#endif
