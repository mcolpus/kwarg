/*******************************************************************

    arg_builder.h

    Description of functions for building and outputting an ancestral
    recombination graph

********************************************************************/

#ifndef ARG_BUILDER_H
#define ARG_BUILDER_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>

typedef enum
{
    GENE_ANY,
    GENE_BEAGLE,
    GENE_FASTA
} Gene_Format;

typedef enum
{
    GENE_BINARY,
    GENE_NUCLEIC,
    GENE_AMINO
} Gene_SeqType;

typedef enum
{
    ARGDOT, /* Output ARG as ARG in DOT format */
    ARGGDL, /* Output ARG as ARG in GDL format */
    ARGGML  /* Output ARG as ARG in GML format */
} ARGOutputFormat;

typedef enum
{
    LABEL_NONE,     /* Do not label nodes */
    LABEL_BOTH,     /* Label nodes with both label and sequence */
    LABEL_SEQUENCE, /* Label nodes only with sequence */
    LABEL_LABEL     /* Label nodes only with label */
} ARGOutputLabels;

typedef enum
{
    UNSET,         /* Node type has not been set yet */
    SAMPLE,        /* Node represents a sampled sequence */
    COALESCENCE,   /* Node represents a coalescence */
    RECOMBINATION, /* Node represents a recombination */
    ROOT           /* Node representing the root */
} NodeType;

typedef struct _Node Node;

typedef struct _Edge
{
    Node *from;
    Node *to;
    std::vector<int> mutations; /* Vector of sites which mutated along this edge */
    std::vector<int> back_mutations;
} Edge;

typedef struct _Node
{
    int id; // This is used to identify the node, particularly useful for printing graph
    NodeType type;
    std::string label; // This is a description of the node. Has no effect on logic
    std::vector<int> mutations; // Most sequences have very few mutations, so easier to store as a vector
    union U
    {
        bool none;
        Edge *one;
        struct
        {
            Edge *prefix;
            Edge *suffix;
            int position;
        } two;

        U()
        {
        }

        ~U()
        {
        }
    } predecessor;
} Node;

typedef struct _ARG
{
    // Unique pointers are used to manage memory
    std::vector<std::unique_ptr<Node>> nodes;
    std::vector<std::unique_ptr<Edge>> edges;
    std::multimap<int, Edge *> mutations_to_edges;
    Node root;
    int number_of_ancestral_nodes = 0;
    int number_of_back_mutations = 0;
    int number_of_recurrent_mutations = 0;
    std::set<int> back_and_recurrent_mutations;

    _ARG()
    {
        root.type = ROOT;
        root.id = -1;
        root.mutations.clear();
        root.label = "Root";
        root.predecessor.none = true;
    }
} ARG;

typedef struct _GENE
{
    std::string label;
    std::vector<int> mutations;
} Gene;

typedef struct _GENEs
{
    std::vector<Gene> genes;
} Genes;

void arg_output(const ARG &arg, const Genes &genes, FILE *fp,
                ARGOutputFormat format, bool annotate_edges, ARGOutputLabels node_labels);
ARG build_arg(Genes genes, FILE *print_progress, int how_verbose);

#endif