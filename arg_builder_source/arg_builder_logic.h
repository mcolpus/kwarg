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
    std::string label;          // This is a description of the node. Has no effect on logic
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
    std::multimap<int, Edge *> mutation_to_edges;          // map to all edges which have the mutation
    std::multimap<int, Edge *> back_mutation_to_edges;     // map to all edges which have the back mutation
    std::multimap<int, Edge *> mutation_to_recombinations; // map to all recombination node out-edges containing the mutation
    // Node root;
    Node *root_ptr;
    std::vector<int> root_mutations;

    int number_of_ancestral_nodes = 0;
    int number_of_back_mutations = 0;
    int number_of_recurrent_mutations = 0;
    int number_of_recombinations = 0;
    std::set<int> recurrent_sites;
    std::set<int> back_mutation_sites;
    std::set<Node *> recombination_nodes;

    std::vector<std::unique_ptr<Edge>> self_edges;      // These are used for certain nodes (roots and recombinations) which are always relevant

    Edge * create_self_loop(Node *node_ptr)
    {
        auto self_edge = std::make_unique<Edge>();
        self_edge->from = node_ptr;
        self_edge->to = node_ptr;

        Edge *edge_ptr = self_edge.get();

        self_edges.push_back(std::move(self_edge));

        return edge_ptr;
    }

    _ARG()
    {
        auto root = std::make_unique<Node>();

        root->type = ROOT;
        root->id = -1;
        root->mutations.clear();
        root->label = "Root";
        root->predecessor.none = true;

        create_self_loop(root.get());

        root_ptr = root.get();
        nodes.push_back(std::move(root));
    }

    _ARG(std::string root_label, std::vector<int> _root_mutations)
    {
        auto root = std::make_unique<Node>();

        root->type = ROOT;
        root->id = -1;
        root->mutations = _root_mutations;
        root->label = root_label;
        root->predecessor.none = true;

        create_self_loop(root.get());

        root_mutations = _root_mutations;
        root_ptr = root.get();

        nodes.push_back(std::move(root));
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

float get_cost(const int rms, const int bms, const int rcs);
float get_cost(const int rms, const int bms);
void arg_output(const ARG &arg, const Genes &genes, FILE *fp,
                ARGOutputFormat format, bool annotate_edges, ARGOutputLabels node_labels);
ARG build_arg(Genes genes, bool root_given, int how_verbose, float cost_rm, float cost_bm, float cost_recomb, int recomb_max, int rm_max, int bm_max);
ARG build_arg_multi_random_runs(int number_of_runs, int run_seed, Genes genes, bool root_given, int how_verbose, float cost_rm, float cost_bm,
                                float cost_recomb, int recomb_max, int rm_max, int bm_max);
ARG build_arg_multi_smart_runs(int number_of_runs, int run_seed, bool tricky_first, const Genes genes, bool root_given, int how_verbose, float cost_rm,
                               float cost_bm, float cost_recomb, int recomb_max, int rm_max, int bm_max);

#endif