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

#include "arg_builder_logic.h"
#include "vector_set_operations.h"

FILE *print_progress;

void print_gene(const Gene &g)
{
    if (g.label != "")
    {
        fprintf(print_progress, g.label.c_str());
        fprintf(print_progress, ": ");
    }

    fprintf(print_progress, "mutations: ");
    for (int i : g.mutations)
    {
        fprintf(print_progress, " %d,", i);
    }
    fprintf(print_progress, "\n");
}

void print_arg(const ARG &arg)
{
    fprintf(print_progress, "--------------printing arg----------------\n");
    fprintf(print_progress, "-----------printing nodes:\n");
    for (auto &node_ptr : arg.nodes)
    {
        if (node_ptr->label != "")
        {
            fprintf(print_progress, node_ptr->label.c_str());
            fprintf(print_progress, ":  ");
        }

        fprintf(print_progress, "mutations: ");
        for (int i : node_ptr->mutations)
        {
            fprintf(print_progress, " %d,", i);
        }
        fprintf(print_progress, "\n");
    }

    fprintf(print_progress, "-----------printing edges:\n");
    for (auto &edge_ptr : arg.edges)
    {
        fprintf(print_progress, "from: %s, to: %s, muts:", edge_ptr->from->label.c_str(), edge_ptr->to->label.c_str());
        for (int i : edge_ptr->mutations)
        {
            fprintf(print_progress, " %d,", i);
        }
        fprintf(print_progress, " back_muts:");
        for (int i : edge_ptr->back_mutations)
        {
            fprintf(print_progress, " %d,", i);
        }
        fprintf(print_progress, "\n");
    }

    fprintf(print_progress, "-----------printing edges from nodes:\n");
    for (auto &node_ptr : arg.nodes)
    {
        if (node_ptr->label != "")
        {
            fprintf(print_progress, node_ptr->label.c_str());
            fprintf(print_progress, ":  ");
        }

        /* Output edges out of this node */
        if (node_ptr->type == SAMPLE || node_ptr->type == COALESCENCE)
        {
            fprintf(print_progress, "%d -> %d,", node_ptr->predecessor.one->from->id, node_ptr->id);
        }
        else if (node_ptr->type == RECOMBINATION)
        {
            fprintf(print_progress, "%d -> %d, and %d -> %d", node_ptr->predecessor.two.prefix->from->id, node_ptr->id,
                    node_ptr->predecessor.two.suffix->from->id, node_ptr->id);
        }

        fprintf(print_progress, "\n");
    }
}

std::string vector_to_string(const std::vector<int> &v, bool make_negative = false)
{
    std::string s = "";
    bool first = true;
    for(const int i : v)
    {
        if(!first)
            s += ", ";
        if (make_negative)
            s += std::to_string(-i);
        else
            s += std::to_string(i);
        first = false;
    }

    return s;
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
        fprintf(fp, "digraph ARG {\n  { rank = same;");
        for (auto const& n : arg.nodes)
        {
            if (n->type == SAMPLE)
                fprintf(fp, " %d;", n->id);
        }
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
            label_printed= true;
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
                fprintf(fp, "  %d -> %d [label=\"P",
                        prefix_id, node_id);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"P",
                        prefix_id, node_id);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"P",
                        prefix_id, node_id);
                break;
            }
            if (annotate_edges)
            {
                auto muts = vector_to_string(arg.nodes[i]->predecessor.two.prefix->mutations);
                auto back_muts = vector_to_string(arg.nodes[i]->predecessor.two.prefix->back_mutations, true);
                fprintf(fp, " %s|%s", muts, back_muts);
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
                fprintf(fp, "  %d -> %d [label=\"S",
                        suffix_id, node_id);
                break;
            case ARGGDL:
                fprintf(fp,
                        "  edge: { sourcename: \"%d\" targetname: \"%d\" label: \"S",
                        suffix_id, node_id);
                break;
            case ARGGML:
                fprintf(fp,
                        "  edge [\n    linestyle \"solid\"\n    source %d\n    target %d\n    label \"S",
                        suffix_id, node_id);
                break;
            }
            if (annotate_edges)
            {
                auto muts = vector_to_string(arg.nodes[i]->predecessor.two.suffix->mutations);
                auto back_muts = vector_to_string(arg.nodes[i]->predecessor.two.suffix->back_mutations, true);
                fprintf(fp, " %s|%s", muts, back_muts);
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

void add_seq_to_arg_rm_only_all_samples(ARG &arg, const Gene &g)
{
    // Build up set of nodes referenced by g's mutations
    std::set<Node *> related_nodes;
    related_nodes.insert(&arg.root);

    std::vector<int> existing_mutations;
    for (int mut : g.mutations)
    {
        auto range = arg.mutations_to_edges.equal_range(mut);
        bool contained = false;
        for (auto i = range.first; i != range.second; ++i)
        {
            related_nodes.insert(i->second->to);
            contained = true;
        }

        if (contained)
            existing_mutations.push_back(mut); // These should be accounted for or else RM
    }

    // TODO: currently inserts via greedy alg to minimize RM's

    int best_score = INT32_MAX;
    Node *best_node = nullptr;
    for (Node *node : related_nodes)
    {
        std::vector<int> recurrent_mutations = (existing_mutations, node->mutations);
        std::vector<int> back_mutations = vector_difference(node->mutations, existing_mutations);

        int score = recurrent_mutations.size() + back_mutations.size();

        if (score < best_score)
        {
            best_score = score;
            best_node = node;
        }
    }

    // Now insert g into a node, child of best node

    auto new_node = std::make_unique<Node>();
    auto new_edge = std::make_unique<Edge>();
    new_node->label = g.label;
    new_node->mutations = g.mutations;
    new_node->type = SAMPLE;
    new_node->predecessor.one = new_edge.get();
    new_edge->from = best_node;
    new_edge->to = new_node.get();

    std::vector<int> new_mutations = vector_difference(g.mutations, best_node->mutations);
    std::vector<int> back_mutations = vector_difference(best_node->mutations, g.mutations);
    for (int m : new_mutations)
        arg.mutations_to_edges.insert({m, new_edge.get()});
    for (int m : back_mutations)
        arg.mutations_to_edges.insert({m, new_edge.get()});

    new_edge->mutations = std::move(new_mutations);
    new_edge->back_mutations = std::move(back_mutations);

    arg.edges.push_back(std::move(new_edge));
    arg.nodes.push_back(std::move(new_node));
}

void add_seq_to_arg_rm_only(ARG &arg, const Gene &g)
{
    // Build up set of edges referenced by g's mutations
    std::set<Edge *> related_edges;

    std::vector<int> existing_mutations;
    for (int mut : g.mutations)
    {
        auto range = arg.mutations_to_edges.equal_range(mut);
        bool contained = false;
        for (auto i = range.first; i != range.second; ++i)
        {
            related_edges.insert(i->second);
            contained = true;
        }

        if (contained)
            existing_mutations.push_back(mut); // These should be accounted for or else RM
    }

    // TODO: currently inserts via greedy alg to minimize RM's

    int best_score = existing_mutations.size(); // This is by just attaching to root
    Edge *best_edge = nullptr;
    for (Edge *edge : related_edges)
    {
        // First find needed mutations relative to target node of edge
        std::vector<int> recurrent_mutations = vector_difference(existing_mutations, edge->to->mutations); // TODO: a symmetric difference is going on
        std::vector<int> back_mutations = vector_difference(edge->to->mutations, existing_mutations);

        // Now consider if some can be removed by splitting edge
        recurrent_mutations = vector_difference(recurrent_mutations, edge->back_mutations);
        back_mutations = vector_difference(back_mutations, edge->mutations);

        int score = recurrent_mutations.size() + back_mutations.size();

        if (score < best_score)
        {
            best_score = score;
            best_edge = edge;
        }
    }

    if (best_edge == nullptr)
    {
        // add to root
        auto new_node = std::make_unique<Node>();
        auto new_edge = std::make_unique<Edge>();

        new_node->id = arg.nodes.size();
        new_node->label = g.label;
        new_node->mutations = g.mutations;
        new_node->type = SAMPLE;
        new_node->predecessor.one = new_edge.get();

        new_edge->from = &arg.root;
        new_edge->to = new_node.get();

        for (int m : g.mutations)
            arg.mutations_to_edges.insert({m, new_edge.get()});

        new_edge->mutations = g.mutations;
        new_edge->back_mutations.clear();

        arg.edges.push_back(std::move(new_edge));
        arg.nodes.push_back(std::move(new_node));

        return;
    }

    // Now insert g into a node, child of best node

    // Need to check if best required splitting edge
    std::vector<int> recurrent_mutations = vector_difference(existing_mutations, best_edge->to->mutations);
    std::vector<int> back_mutations = vector_difference(best_edge->to->mutations, existing_mutations);

    auto [removable_rms, required_rms] = vector_split(recurrent_mutations, best_edge->back_mutations);
    auto [removable_back_muts, required_back_muts] = vector_split(back_mutations, best_edge->mutations);
    auto new_mutations = vector_difference(g.mutations, existing_mutations);

    if (removable_rms.size() > 0 || removable_back_muts.size() > 0)
    {
        // Should split, but is splitting enough?
        // 1. Split edge first
        auto split_node = std::make_unique<Node>(); // Between the top and bottom of best_edge
        auto split_edge = std::make_unique<Edge>(); // From top to split_node

        split_edge->from = best_edge->from;
        split_edge->to = best_edge->from = split_node.get();

        split_edge->mutations = vector_difference(best_edge->mutations, removable_back_muts);
        split_edge->back_mutations = vector_difference(best_edge->back_mutations, removable_rms);
        best_edge->mutations = removable_back_muts;
        best_edge->back_mutations = removable_rms;

        split_node->id = arg.nodes.size();
        split_node->label = "ancestral" + std::to_string(arg.number_of_ancestral_nodes);
        split_node->type = COALESCENCE;
        split_node->predecessor.one = split_edge.get();
        split_node->mutations = vector_difference(vector_union(split_edge->from->mutations, split_edge->mutations), split_edge->back_mutations);

        for (int m : split_edge->mutations)
        {
            arg.mutations_to_edges.insert({m, split_edge.get()});
            remove_edge_from_multimap(arg.mutations_to_edges, m, best_edge);
        }
        for (int m : split_edge->back_mutations)
        {
            arg.mutations_to_edges.insert({m, split_edge.get()});
            remove_edge_from_multimap(arg.mutations_to_edges, m, best_edge);
        }

        // Now check if the split node is the node we wish to insert
        if (split_node->mutations == g.mutations)
        {
            // Then relabel split
            split_node->label = g.label;
            split_node->type = SAMPLE;

            arg.edges.push_back(std::move(split_edge));
            arg.nodes.push_back(std::move(split_node));
        }
        else
        {
            // Add new node as child of split_node
            auto new_node = std::make_unique<Node>();
            auto new_edge = std::make_unique<Edge>(); // From split_node to new_node

            new_edge->from = split_node.get();
            new_edge->to = new_node.get();

            new_edge->mutations = vector_union(required_rms, vector_difference(g.mutations, existing_mutations));
            new_edge->back_mutations = required_back_muts;

            new_node->id = arg.nodes.size() + 1; // split_node has id before
            new_node->label = g.label;
            new_node->type = SAMPLE;
            new_node->predecessor.one = new_edge.get();
            new_node->mutations = g.mutations;

            for (int m : new_edge->mutations)
                arg.mutations_to_edges.insert({m, new_edge.get()});
            for (int m : new_edge->back_mutations)
                arg.mutations_to_edges.insert({m, new_edge.get()});

            arg.edges.push_back(std::move(new_edge));
            arg.edges.push_back(std::move(split_edge));
            arg.nodes.push_back(std::move(new_node));
            arg.nodes.push_back(std::move(split_node));
            arg.number_of_ancestral_nodes += 1;
        }
    }
    else
    {
        // Add new_node as child of target of best_edge
        auto new_node = std::make_unique<Node>();
        auto new_edge = std::make_unique<Edge>();

        new_node->id = arg.nodes.size();
        new_node->label = g.label;
        new_node->mutations = g.mutations;
        new_node->type = SAMPLE;
        new_node->predecessor.one = new_edge.get();

        new_edge->from = best_edge->to;
        new_edge->to = new_node.get();

        std::vector<int> new_mutations = vector_difference(g.mutations, best_edge->to->mutations);
        std::vector<int> back_mutations = vector_difference(best_edge->to->mutations, g.mutations);
        for (int m : new_mutations)
            arg.mutations_to_edges.insert({m, new_edge.get()});
        for (int m : back_mutations)
            arg.mutations_to_edges.insert({m, new_edge.get()});

        new_edge->mutations = std::move(new_mutations);
        new_edge->back_mutations = std::move(back_mutations);

        arg.edges.push_back(std::move(new_edge));
        arg.nodes.push_back(std::move(new_node));
    }
}

ARG build_arg(Genes genes, FILE *in_print_progress, bool print_steps = false)
{
    print_progress = in_print_progress;
    fprintf(print_progress, "starting build\n");
    ARG arg;
    int step = 0;
    for (Gene &g : genes.genes)
    {
        add_seq_to_arg_rm_only(arg, g);
        step += 1;
        if (print_steps)
        {
            std::string filename = "dotty";
            filename += std::to_string(step) + ".dot";
            auto fp = fopen(filename.c_str(), "w");
            arg_output(arg, genes, fp, ARGDOT, true, LABEL_BOTH);
            fclose(fp);
        }
    }

    print_arg(arg);

    fprintf(print_progress, "ARG built\n");

    return arg;
}
