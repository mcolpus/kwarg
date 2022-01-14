/***************************************************************************
 *
 *    arg_builder.cpp
 *
 *    Implementation of logic for building ARGs by joining one sequence at
 *    a time.
 *
 ****************************************************************************/

#include <vector>
#include <string>
#include <map>
#include <set>
#include <memory>
#include <algorithm>
#include <iostream>

#include "arg_builder_logic.h"
#include "vector_set_operations.h"

int _how_verbose = 0;
float _cost_recurrent_mutation = 1.0;
float _cost_back_mutation = 1.0;
float _cost_recombination = 1.0;

std::string vector_to_string(const std::vector<int> &v, bool make_negative = false)
{
    std::string s = "";
    bool first = true;
    for (const int i : v)
    {
        if (!first)
            s += ", ";
        if (make_negative)
            s += std::to_string(-i);
        else
            s += std::to_string(i);
        first = false;
    }

    return s;
}

void print_gene(const Gene &g)
{
    if (g.label != "")
    {
        std::cout << g.label << ": ";
    }

    std::cout << "mutations: ";
    for (int i : g.mutations)
    {
        std::cout << i << ", ";
    }
    std::cout << "\n";
}

void print_arg(const ARG &arg)
{
    std::cout << "\n\n--------------printing arg----------------\n\n";
    std::cout << "-----------printing nodes:\n";
    for (auto &node_ptr : arg.nodes)
    {
        if (node_ptr->label != "")
        {
            std::cout << node_ptr->label << ": ";
        }

        std::cout << "mutations: ";
        for (int i : node_ptr->mutations)
        {
            std::cout << i << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "-----------printing edges:\n";
    for (auto &edge_ptr : arg.edges)
    {
        std::cout << edge_ptr->from->label << " -> " << edge_ptr->to->label << " muts: ";
        for (int i : edge_ptr->mutations)
        {
            std::cout << i << ", ";
        }
        std::cout << " back_muts:";
        for (int i : edge_ptr->back_mutations)
        {
            std::cout << i << ", ";
        }
        std::cout << "\n";
    }

    std::cout << "-----------printing edges from nodes:\n";
    for (auto &node_ptr : arg.nodes)
    {
        if (node_ptr->label != "")
        {
            std::cout << node_ptr->label << ": ";
        }

        /* Output edges out of this node */
        if (node_ptr->type == SAMPLE || node_ptr->type == COALESCENCE)
        {
            std::cout << node_ptr->predecessor.one->from->id << " -> " << node_ptr->id;
        }
        else if (node_ptr->type == RECOMBINATION)
        {
            std::cout << node_ptr->predecessor.two.prefix->from->id << " -> " << node_ptr->id << " <- " << node_ptr->predecessor.two.suffix->from->id;
        }
        std::cout << "\n";
    }
}

float get_cost(const int rms, const int bms, const int rcs)
{
    if (_cost_recurrent_mutation < 0 && rms > 0)
        return -1.0;
    else if (_cost_back_mutation < 0 && bms > 0)
        return -1.0;
    else if (_cost_recombination < 0 && rcs > 0)
        return -1.0;

    return static_cast<float>(rms) * _cost_recurrent_mutation +
           static_cast<float>(bms) * _cost_back_mutation +
           static_cast<float>(rcs) * _cost_recombination;
}

float get_cost(const int rms, const int bms)
{
    return get_cost(rms, bms, 0);
}

void arg_output(const ARG &arg, const Genes &genes, FILE *fp,
                ARGOutputFormat format, bool annotate_edges, ARGOutputLabels node_labels)
{
    int intervals;

    if (fp == NULL)
        fp = stdout;

    /* Output prelude */
    switch (format)
    {
    case ARGDOT:
        fprintf(fp, "digraph ARG {\n");
        break;
    case ARGGDL:
        fprintf(fp, "graph: {\n");
        break;
    case ARGGML:
        fprintf(fp, "graph [\n  directed 1\n");
        break;
    }

    // Output the root node
    switch (format)
    {
    case ARGDOT:
        fprintf(fp, "  -1 [label=\"-1 root");
        break;
    case ARGGDL:
        fprintf(fp, "  node: { title: \"-1\" label: \"-1 root");
        break;
    case ARGGML:
        fprintf(fp, "  node [\n    id -1\n    label \"-1 root");
        break;
    }
    fprintf(fp, "\"");

    switch (format)
    {
    case ARGDOT:
        fprintf(fp, ",color=black];\n");
        break;
    case ARGGDL:
        fprintf(fp, " shape: circle bordercolor: black }\n");
        break;
    case ARGGML:
        fprintf(fp, "\n    graphics [\n      outline \"black");
        fprintf(fp, "\"\n    ]\n    vgj [\n      type \"Oval\"\n      labelPosition \"in\"\n    ]\n  ]\n");
        break;
    }

    /* Output nodes and their edges */
    for (int i = 0; i < arg.nodes.size(); i++)
    {
        int node_id = arg.nodes[i]->id;

        /* Output node i */
        bool label_printed = false;
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, "  %d [label=\"", node_id);
            break;
        case ARGGDL:
            fprintf(fp, "  node: { title: \"%d\" label: \"", node_id);
            break;
        case ARGGML:
            fprintf(fp, "  node [\n    id %d\n    label \"", node_id);
            break;
        }

        /* Generate node label */
        if (arg.nodes[i]->label != "" && (node_labels == LABEL_BOTH || node_labels == LABEL_LABEL))
        {
            fprintf(fp, "%s", arg.nodes[i]->label.c_str());
            label_printed = true;
        }
        if (node_labels == LABEL_BOTH || node_labels == LABEL_SEQUENCE)
        {
            if (label_printed)
                fprintf(fp, "; %s", vector_to_string(arg.nodes[i]->mutations).c_str());
            else
                fprintf(fp, "%s", vector_to_string(arg.nodes[i]->mutations).c_str());
            label_printed = true;
        }
        if (!label_printed)
        {
            /* No node label printed */
            if ((arg.nodes[i]->type == RECOMBINATION) && (arg.nodes[i]->label != ""))
                /* If a recombination node is still unlabeled, label it (label
                 * should be recombination point) if possible.
                 */
                fprintf(fp, "%s\"", arg.nodes[i]->label.c_str());
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
        switch (arg.nodes[i]->type)
        {
        case SAMPLE:
            fprintf(fp, "red");
            break;
        case COALESCENCE:
            fprintf(fp, "green");
            break;
        case RECOMBINATION:
            fprintf(fp, "blue");
            break;
        case ROOT:
            fprintf(fp, "black");
            break;
        case UNSET:
            fprintf(fp, "purple");
            break;
        }
        /* End node description */
        switch (format)
        {
        case ARGDOT:
            fprintf(fp, "];\n");
            break;
        case ARGGDL:
            if (i < genes.genes.size())
                fprintf(fp, " vertical_order: maxlevel");
            fprintf(fp, " }\n");
            break;
        case ARGGML:
            fprintf(fp,
                    "\"\n    ]\n    vgj [\n      type \"Oval\"\n      labelPosition \"in\"\n    ]\n  ]\n");
            break;
        }

        /* Output edges out of this node */
        switch (arg.nodes[i]->type)
        {
        case SAMPLE:
        case COALESCENCE:
        {
            /* One outgoing edge */
            int from_id = arg.nodes[i]->predecessor.one->from->id;
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d", from_id, node_id);
                break;
            case ARGGDL:
                fprintf(fp, "  edge: { sourcename: \"%d\" targetname: \"%d\"",
                        from_id, node_id);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d",
                        from_id, node_id);
                break;
            }
            /* Output edge label */
            if (annotate_edges)
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
                auto muts = vector_to_string(arg.nodes[i]->predecessor.one->mutations);
                auto back_muts = vector_to_string(arg.nodes[i]->predecessor.one->back_mutations, true);
                fprintf(fp, "%s|%s", muts.c_str(), back_muts.c_str());
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
        }
        break;
        case RECOMBINATION:
        {
            /* Two outgoing edges */
            /* Prefix edge */
            int prefix_id = arg.nodes[i]->predecessor.two.prefix->from->id;
            int suffix_id = arg.nodes[i]->predecessor.two.suffix->from->id;
            int crossover_position = arg.nodes[i]->predecessor.two.position;
            switch (format)
            {
            case ARGDOT:
                fprintf(fp, "  %d -> %d [label=\"P %d",
                        prefix_id, node_id, crossover_position);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"P %d",
                        prefix_id, node_id, crossover_position);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"P %d",
                        prefix_id, node_id, crossover_position);
                break;
            }
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
                fprintf(fp, "  %d -> %d [label=\"S %d",
                        suffix_id, node_id, crossover_position);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"S %d",
                        suffix_id, node_id, crossover_position);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"S %d",
                        suffix_id, node_id, crossover_position);
                break;
            }
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
        }
        break;
        case ROOT:
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

void remove_edge_from_multimap(std::multimap<int, Edge *> &map, int mut, Edge *edge)
{
    auto iterpair = map.equal_range(mut);

    for (auto it = iterpair.first; it != iterpair.second; ++it)
    {
        if (it->second == edge)
        {
            map.erase(it);
            break;
        }
    }
}

std::tuple<std::set<Edge *>, std::vector<int>> find_relevant_edges(ARG &arg, const Gene &g)
{
    std::set<Edge *> relevant_edges;

    std::vector<int> existing_mutations;
    for (int mut : g.mutations)
    {
        // Any edge containing a mutation which is in g may be relevant
        auto range = arg.mutation_to_edges.equal_range(mut);
        bool contained = false;
        for (auto i = range.first; i != range.second; ++i)
        {
            relevant_edges.insert(i->second);
            contained = true;
        }
        if (contained)
            existing_mutations.push_back(mut); // These should be accounted for or else RM

        // Any recombination node which shares mutations with g may be relevant
        range = arg.mutation_to_recombinations.equal_range(mut);
        for (auto i = range.first; i != range.second; ++i)
        {
            relevant_edges.insert(i->second);
        }
    }

    // Any edge with a back mutation for a site not in g may be relevant.
    std::vector<int> relevant_bms;
    std::set_difference(arg.back_mutation_sites.begin(), arg.back_mutation_sites.end(),
                        g.mutations.begin(), g.mutations.end(),
                        std::inserter(relevant_bms, relevant_bms.begin()));
    for (int back_mut : relevant_bms)
    {
        auto range = arg.back_mutation_to_edges.equal_range(back_mut);
        for (auto i = range.first; i != range.second; ++i)
        {
            relevant_edges.insert(i->second);
        }
    }

    return std::make_tuple(relevant_edges, existing_mutations);
}

/* Splits the edge into two and returns the new node. removable mutations are put below, rest are above */
Node *split_edge(ARG &arg, Edge *edge, std::vector<int> removable_muts, std::vector<int> removable_back_muts)
{
    auto split_node = std::make_unique<Node>(); // Between the top and bottom of best_edge
    auto split_edge = std::make_unique<Edge>(); // From top to split_node

    split_edge->from = edge->from;
    split_edge->to = edge->from = split_node.get();

    split_edge->mutations = vector_difference(edge->mutations, removable_back_muts);
    split_edge->back_mutations = vector_difference(edge->back_mutations, removable_muts);
    edge->mutations = removable_back_muts;
    edge->back_mutations = removable_muts;

    split_node->id = arg.nodes.size();
    split_node->label = "ancestral" + std::to_string(arg.number_of_ancestral_nodes);
    split_node->type = COALESCENCE;
    split_node->predecessor.one = split_edge.get();
    split_node->mutations = vector_difference(vector_union(split_edge->from->mutations, split_edge->mutations), split_edge->back_mutations);

    for (int m : split_edge->mutations)
    {
        arg.mutation_to_edges.insert({m, split_edge.get()});
        remove_edge_from_multimap(arg.mutation_to_edges, m, edge);
    }
    for (int m : split_edge->back_mutations)
    {
        arg.back_mutation_to_edges.insert({m, split_edge.get()});
        remove_edge_from_multimap(arg.back_mutation_to_edges, m, edge);
    }

    Node *new_node = split_node.get();

    arg.edges.push_back(std::move(split_edge));
    arg.nodes.push_back(std::move(split_node));
    arg.number_of_ancestral_nodes += 1;

    return new_node;
}

Node *insert_seq_as_direct_child(ARG &arg, const Gene &g, Node *parent, const std::vector<int> &existing_mutations)
{
    auto new_node = std::make_unique<Node>();
    auto new_edge = std::make_unique<Edge>();

    new_node->id = arg.nodes.size();
    new_node->label = g.label;
    new_node->mutations = g.mutations;
    new_node->type = SAMPLE;
    new_node->predecessor.one = new_edge.get();

    new_edge->from = parent;
    new_edge->to = new_node.get();

    new_edge->mutations = vector_difference(g.mutations, parent->mutations);
    new_edge->back_mutations = vector_difference(parent->mutations, g.mutations);
    auto recurrent_muts = vector_intersect(new_edge->mutations, existing_mutations);

    for (int m : new_edge->mutations)
        arg.mutation_to_edges.insert({m, new_edge.get()});
    for (int m : new_edge->back_mutations)
        arg.back_mutation_to_edges.insert({m, new_edge.get()});

    // All the existing mutations are now count towards recurrent mutations
    arg.number_of_recurrent_mutations += recurrent_muts.size();
    arg.number_of_back_mutations += new_edge->back_mutations.size();
    for (auto m : recurrent_muts)
    {
        arg.recurrent_sites.insert(m);
    }
    for (auto m : new_edge->back_mutations)
    {
        arg.back_mutation_sites.insert(m);
    }

    // If the parent node is a recombination node then we want to update arg.mutation_to_recombinations;
    if (parent->type == RECOMBINATION)
    {
        for (int m : parent->mutations)
        {
            arg.mutation_to_recombinations.insert({m, new_edge.get()});
        }
    }

    Node *return_value = new_node.get();

    arg.edges.push_back(std::move(new_edge));
    arg.nodes.push_back(std::move(new_node));

    return return_value;
}

Node *recombine_nodes(ARG &arg, const int pos, Node *prefix, Node *suffix)
{
    auto recomb_node = std::make_unique<Node>();
    auto prefix_edge = std::make_unique<Edge>();
    auto suffix_edge = std::make_unique<Edge>();

    prefix_edge->from = prefix;
    prefix_edge->to = recomb_node.get();
    suffix_edge->from = suffix;
    suffix_edge->to = recomb_node.get();

    recomb_node->id = arg.nodes.size();
    recomb_node->label = "recomb" + std::to_string(arg.number_of_ancestral_nodes);
    recomb_node->mutations = vector_union(vector_values_below(prefix->mutations, pos),
                                          vector_values_above(suffix->mutations, pos));
    recomb_node->type = RECOMBINATION;
    recomb_node->predecessor.two.position = pos;
    recomb_node->predecessor.two.prefix = prefix_edge.get();
    recomb_node->predecessor.two.suffix = suffix_edge.get();

    Node *return_value = recomb_node.get();
    arg.recombination_nodes.insert(return_value);

    arg.nodes.push_back(std::move(recomb_node));
    arg.edges.push_back(std::move(prefix_edge));
    arg.edges.push_back(std::move(suffix_edge));
    arg.number_of_ancestral_nodes += 1;
    arg.number_of_recombinations += 1;

    return return_value;
}

std::tuple<Edge *, int, int> find_best_single_parent_location(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    // set best score initially to correspond to attaching to root
    int best_rms = static_cast<float>(existing_mutations.size());
    int best_bms = 0;
    float best_score = get_cost(best_rms, best_bms); // Could be negative if recurrent mutations aren't allowed
    Edge *best_edge = nullptr;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        back_mutations = vector_difference(back_mutations, edge->mutations);

        float score = get_cost(recurrent_mutations.size(), back_mutations.size());

        if (score >= 0 && (score < best_score || best_score < 0))
        {
            best_score = score;
            best_rms = recurrent_mutations.size();
            best_bms = back_mutations.size();
            best_edge = edge;
        }
    }

    return std::make_tuple(best_edge, best_rms, best_bms);
}

void insert_seq_at_edge(ARG &arg, const Gene &g, Edge *best_edge, const std::vector<int> &existing_mutations)
{
    if (best_edge == nullptr)
    {
        // Insert at root
        insert_seq_as_direct_child(arg, g, &arg.root, existing_mutations);
        return;
    }

    // Need to check if best required splitting edge
    std::vector<int> recurrent_mutations = vector_difference(existing_mutations, best_edge->to->mutations);
    std::vector<int> back_mutations = vector_difference(best_edge->to->mutations, existing_mutations);

    auto [removable_rms, required_rms] = vector_split(recurrent_mutations, best_edge->back_mutations);
    auto [removable_back_muts, required_back_muts] = vector_split(back_mutations, best_edge->mutations);
    auto new_mutations = vector_difference(g.mutations, existing_mutations);

    if (removable_rms.size() > 0 || removable_back_muts.size() > 0)
    {
        Node *split_node = split_edge(arg, best_edge, removable_rms, removable_back_muts);

        // Now check if the split node is the node we wish to insert
        if (split_node->mutations == g.mutations)
        {
            // Then relabel split
            split_node->label = g.label;
            split_node->type = SAMPLE;
            arg.number_of_ancestral_nodes -= 1;
        }
        else
        {
            insert_seq_as_direct_child(arg, g, split_node, existing_mutations);
        }
    }
    else
    {
        insert_seq_as_direct_child(arg, g, best_edge->to, existing_mutations);
    }
}

void replace_cost_if_better(std::map<int, std::tuple<Edge *, int, int>> &costs, int pos, Edge *edge, int rms, int bms)
{
    // add entry to prefix_costs
    auto search = costs.find(pos);
    if (search != costs.end())
    {
        auto [edge, current_rms, current_bms] = search->second;
        float new_cost = get_cost(rms, bms);
        if (new_cost >= 0 && new_cost < get_cost(current_rms, current_bms))
            costs.insert_or_assign(pos, std::make_tuple(edge, rms, bms));
    }
    else
    {
        costs.insert({pos, std::make_tuple(edge, rms, bms)});
    }
}

std::tuple<int, Edge *, Edge *, int, int> find_best_recomb_location(ARG &arg, const Gene &g, const std::set<Edge *> &related_edges, const std::vector<int> &existing_mutations)
{
    /**
     * This function will return a list of pairs of edges and associated metrics.
     * Each item in list represents a potential recombination join location.
     * Only considers single cross-over recombination.
     */

    /**
     * For each edge we want to calculate the cost of using it as a prefix or suffix.
     * position 0 refers to just before site 0, 1 just before 1 etc.
     * We don't bother considering the position after the end.
     *
     * prefix_costs is a map from cross-over position to a (edge, rms, bms) pair.
     * suffix_costs will be similar and even have the same keys.
     */
    std::map<int, std::tuple<Edge *, int, int>> prefix_costs;
    std::map<int, std::tuple<Edge *, int, int>> suffix_costs;

    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        auto req_recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        auto req_back_mutations = vector_difference(back_mutations, edge->mutations);

        int i = 0;
        int j = 0;
        int pos = 0;
        int rms = req_recurrent_mutations.size();
        int bms = req_back_mutations.size();

        // set 0 suffix cost and prefix cost
        replace_cost_if_better(suffix_costs, 0, edge, rms, bms);

        while (i < rms && j < bms)
        {
            bool inc_i;
            if (req_recurrent_mutations[i] <= req_back_mutations[j])
            {
                pos = req_recurrent_mutations[i];
                inc_i = true;
            }
            else
            {
                pos = req_back_mutations[j];
                inc_i = false;
            }

            // update prefix and suffix costs
            replace_cost_if_better(prefix_costs, pos, edge, i, j);
            replace_cost_if_better(suffix_costs, pos, edge, rms - i, bms - j);

            if (inc_i)
                i += 1;
            else
                j += 1;
        }

        while (i < rms)
        {
            // So j must equal bms
            pos = req_recurrent_mutations[i];

            // update prefix and suffix costs
            replace_cost_if_better(prefix_costs, pos, edge, i, bms);
            replace_cost_if_better(suffix_costs, pos, edge, rms - i, 0);

            i += 1;
        }
        while (j < bms)
        {
            // So j must equal bms
            pos = req_back_mutations[j];

            // update prefix and suffix costs
            replace_cost_if_better(prefix_costs, pos, edge, rms, j);
            replace_cost_if_better(suffix_costs, pos, edge, 0, bms - j);

            j += 1;
        }

        replace_cost_if_better(suffix_costs, pos + 1, edge, 0, 0);
    }

    // Now we can go through the costs (which should be sorted by key)
    // First turn maps into vectors so that we can iterate through both nicely
    std::vector<std::tuple<int, Edge *, int, int>> prefix_costs_vec;
    std::vector<std::tuple<int, Edge *, int, int>> suffix_costs_vec;
    for (const auto &[site, value] : prefix_costs)
    {
        auto [edge, rms, bms] = value;
        prefix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }
    for (const auto &[site, value] : suffix_costs)
    {
        auto [edge, rms, bms] = value;
        suffix_costs_vec.push_back(std::make_tuple(site, edge, rms, bms));
    }

    // Now go through backwards to create optimal prefix costs, and forwards for optimal suffix costs
    std::vector<std::tuple<int, Edge *, int, int>> optimal_prefix_costs;
    std::vector<std::tuple<int, Edge *, int, int>> optimal_suffix_costs;
    float best_cost = -1.0;
    for (int i = prefix_costs_vec.size() - 1; i >= 0; i--)
    {
        auto [site, edge, rms, bms] = prefix_costs_vec[i];
        float cost = get_cost(rms, bms);
        if (cost > 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_prefix_costs.push_back(prefix_costs_vec[i]);
            best_cost = cost;
        }
    }
    std::reverse(optimal_prefix_costs.begin(), optimal_prefix_costs.end());
    best_cost = -1.0;
    for (int i = 0; i < suffix_costs_vec.size(); i++)
    {
        auto [site, edge, rms, bms] = suffix_costs_vec[i];
        float cost = get_cost(rms, bms);
        if (cost > 0 && (cost < best_cost || best_cost < 0))
        {
            optimal_suffix_costs.push_back(suffix_costs_vec[i]);
            best_cost = cost;
        }
    }

    if (optimal_prefix_costs.size() == 0 or optimal_suffix_costs.size() == 0)
    {
        return std::make_tuple(-1, nullptr, nullptr, 0, 0);
    }

    int p_index = 0, s_index = 0;
    auto [p_site, p_edge, p_rms, p_bms] = optimal_prefix_costs[p_index];
    auto [s_site, s_edge, s_rms, s_bms] = optimal_suffix_costs[s_index];

    int best_pos = -1;
    best_cost = -1.0;
    int best_rms = -1;
    int best_bms = -1;
    Edge *best_prefix = nullptr;
    Edge *best_suffix = nullptr;

    while (p_index < optimal_prefix_costs.size())
    {
        // Now need to find the highest s_index such that
        while (s_index < optimal_suffix_costs.size() - 1)
        {
            int next_site = std::get<0>(optimal_suffix_costs[s_index + 1]);
            if (next_site <= p_site)
            {
                s_index += 1;
                std::tie(s_site, s_edge, s_rms, s_bms) = optimal_suffix_costs[s_index];
            }
            else
            {
                break;
            }
        }

        float current_cost = get_cost(p_rms, p_bms) + get_cost(s_rms, s_bms);
        if (current_cost >= 0 && (current_cost < best_cost || best_cost < 0))
        {
            best_pos = p_site;
            best_cost = current_cost;
            best_rms = p_rms + s_rms;
            best_bms = p_bms + s_bms;
            best_prefix = p_edge;
            best_suffix = s_edge;
        }

        p_index++;
        std::tie(p_site, p_edge, p_rms, p_bms) = optimal_prefix_costs[p_index];
    }

    if (best_prefix == nullptr)
    {
        if (_how_verbose >= 2)
            std::cout << "Could not find recombination location for g: " << g.label << "\n";
        return std::make_tuple(-1, nullptr, nullptr, -1, -1);
    }
    else if (best_suffix == best_prefix)
    {
        if (_how_verbose >= 2)
            std::cout << "Recombination just used one sequence for g: " << g.label << "\n";
        return std::make_tuple(-1, nullptr, nullptr, -1, -1);
    }
    else
    {
        if (_how_verbose >= 2)
            std::cout << "Best recomb is at pos: " << best_pos << ", cost: " << best_cost << ", prefix: " << best_prefix->to->label << ", suffix: " << best_suffix->to->label << "\n";
        return std::make_tuple(best_pos, best_prefix, best_suffix, best_rms, best_bms);
    }
}

void insert_seq_at_recomb_location(ARG &arg, const Gene &g, const int pos, Edge *prefix, Edge *suffix, const std::vector<int> &existing_mutations)
{
    auto p_existing_mutations = vector_values_below(existing_mutations, pos);
    auto s_existing_mutations = vector_values_above(existing_mutations, pos);

    // Check if we should split prefix
    std::vector<int> prefix_rms = vector_difference(p_existing_mutations, prefix->to->mutations);
    std::vector<int> prefix_bms = vector_difference(vector_values_below(prefix->to->mutations, pos),
                                                    p_existing_mutations);

    auto [prefix_removable_rms, prefix_required_rms] = vector_split(prefix_rms, prefix->back_mutations);
    auto [prefix_removable_bms, prefix_required_bms] = vector_split(prefix_bms, prefix->mutations);

    Node *prefix_node; // This will either be the target of the prefix edge, or a split node
    if (prefix_removable_rms.size() > 0 || prefix_removable_bms.size() > 0)
    {
        prefix_node = split_edge(arg, prefix, prefix_removable_rms, prefix_removable_bms);
    }
    else
    {
        prefix_node = prefix->to;
    }

    // Now the same for suffix
    std::vector<int> suffix_rms = vector_difference(s_existing_mutations, suffix->to->mutations);
    std::vector<int> suffix_bms = vector_difference(vector_values_above(suffix->to->mutations, pos),
                                                    s_existing_mutations);

    auto [suffix_removable_rms, suffix_required_rms] = vector_split(suffix_rms, suffix->back_mutations);
    auto [suffix_removable_bms, suffix_required_bms] = vector_split(suffix_bms, suffix->mutations);

    Node *suffix_node; // This will either be the target of the suffix edge, or a split node
    if (suffix_removable_rms.size() > 0 || suffix_removable_bms.size() > 0)
    {
        suffix_node = split_edge(arg, suffix, suffix_removable_rms, suffix_removable_bms);
    }
    else
    {
        suffix_node = suffix->to;
    }

    Node *recomb_node = recombine_nodes(arg, pos, prefix_node, suffix_node);
    insert_seq_as_direct_child(arg, g, recomb_node, existing_mutations);
}

void add_seq_to_arg(ARG &arg, const Gene &g)
{
    auto [relevant_edges, existing_mutations] = find_relevant_edges(arg, g);

    auto [edge, rms, bms] = find_best_single_parent_location(arg, g, relevant_edges, existing_mutations);
    float sp_cost = get_cost(rms, bms);

    if (_cost_recombination >= 0)
    {
        auto [crossover_pos, prefix_edge, suffix_edge, r_rms, r_bms] = find_best_recomb_location(arg, g, relevant_edges, existing_mutations);
        // crossover_pos = -1 means failure, 0 means all one sequence

        float r_cost = get_cost(r_rms, r_bms, 1);

        if (crossover_pos > 0 && r_cost >= 0 && (sp_cost < 0 || r_cost < sp_cost))
        {
            // recombination is valid and better than having single parent
            insert_seq_at_recomb_location(arg, g, crossover_pos, prefix_edge, suffix_edge, existing_mutations);
            return;
        }
    }

    if (sp_cost >= 0)
    {
        insert_seq_at_edge(arg, g, edge, existing_mutations);
    }
    else
    {
        std::cerr << "The gene: " << g.label << " could not be added to the ARG so will be ignored";
    }
}

ARG build_arg(Genes genes, int how_verbose, float cost_rm, float cost_bm, float cost_recomb)
{
    _how_verbose = how_verbose;
    _cost_recurrent_mutation = cost_rm;
    _cost_back_mutation = cost_bm;
    _cost_recombination = cost_recomb;

    if (_how_verbose >= 1)
        std::cout << "starting build\n";

    ARG arg;
    int step = 0;
    for (Gene &g : genes.genes)
    {
        add_seq_to_arg(arg, g);

        step += 1;
        if (_how_verbose >= 3)
        {
            std::string filename = "dotty";
            filename += std::to_string(step) + ".dot";
            auto fp = fopen(filename.c_str(), "w");
            arg_output(arg, genes, fp, ARGDOT, true, LABEL_BOTH);
            fclose(fp);
        }
    }

    if (_how_verbose >= 2)
        print_arg(arg);

    if (_how_verbose >= 1)
    {
        std::cout << "ARG built\n";
        std::cout << "recurrent mutations were: " << vector_to_string(set_to_vector(arg.recurrent_sites)) << "\n";
        std::cout << "back mutations mutations were: " << vector_to_string(set_to_vector(arg.back_mutation_sites)) << "\n";
        std::cout << "ARG made using " << arg.number_of_recurrent_mutations << " recurrent mutations, " << arg.number_of_back_mutations; 
        std::cout << " back mutations, and " << arg.number_of_recombinations << " recombinations\n";
    }

    return arg;
}
