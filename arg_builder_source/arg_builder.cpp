#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <errno.h>

#include <iostream>
#include <fstream>
#include <map>
#include <array>
#include <string>
#include <random>

#include "arg_builder_logic.h"
#include "vector_set_operations.h"

/* Output s to fp with line length l and indentation i */
void pretty_print(FILE *fp, char *s, int l, int i)
{
    int j;
    char *last = s + strlen(s) - 1;

    if (fp == NULL)
        fp = stdout;

    /* Remove initial stretches of white space */
    while (isspace(*s))
        s++;
    if (s > last)
        /* s contains only white space */
        return;
    while (isspace(*last))
        last--;

    while (s <= last)
    {
        /* Output indentation */
        for (j = 0; j < i; j++)
            fputc(' ', fp);

        /* Is there a new line within reach? */
        for (j = 0; (j <= l - i - 2) && (s + j <= last); j++)
            if (s[j] == '\n')
            {
                fwrite(s, sizeof(char), j + 1, fp);
                s += j + 1;
                break;
            }
        if ((j <= l - i - 2) && (s + j <= last))
            /* New line encountered - continue with next line */
            continue;

        if (s + l - i <= last)
        {
            /* Find good place for next line break */
            for (j = l - i - 2; j > 0; j--)
            {
                if ((s[j] == '-') && (s[j - 1] != ' '))
                {
                    /* Break after hyphen */
                    j++;
                    break;
                }
                if ((s[j] == ' ') && (s[j + 1] != '-'))
                {
                    /* Break at space, unless it borders a dash */
                    while ((j > 0) && (s[j] == ' '))
                        j--;
                    if ((j > 0) && (s[j - 1] != '-'))
                        break;
                }
            }
            if (j == 0)
            {
                /* No good line break found - look for acceptable line break */
                for (j = l - i - 2; (j > 0) && (s[j] != ' '); j--)
                    ;
                while (s[j] == ' ')
                    j--;
                if (j == 0)
                    j = l - i - 2;
            }
        }
        else
            /* Rest of text fits on one line */
            j = last - s;
        fwrite(s, sizeof(char), j + 1, fp);
        fputc('\n', fp);
        s += j + 1;
        /* Remove initial stretch of white space */
        while (isspace(*s))
            s++;
    }
}

/* Print an option description to fp with line length l and subsequent
 * line indentation i (length of option if i negative).
 */
void print_option(FILE *fp, char *option, char *description, int l, int i)
{
    int n, m, j;

    if (fp == NULL)
        fp = stdout;

    /* Remove initial stretches of white space */
    while (isspace(*option))
        option++;
    n = strlen(option);
    if (n > 0)
        while (isspace(option[n - 1]))
            n--;
    while (isspace(*description))
        description++;
    m = strlen(description);
    if (m > 0)
        while (isspace(description[m - 1]))
            m--;

    /* Output option */
    fputc(' ', fp);
    fwrite(option, sizeof(char), n, fp);
    fputc(' ', fp);

    if (n + 2 + m <= l)
    {
        /* Whole description fits on first line - output it */
        fwrite(description, sizeof(char), m, fp);
        fputc('\n', fp);
    }
    else
    {
        /* Check whether there is a new line within reach */
        for (j = 0; j <= l - n - 2; j++)
            if (description[j] == '\n')
            {
                fwrite(description, sizeof(char), j + 1, fp);
                pretty_print(fp, description + j + 1, l, (i < 0 ? n + 2 : i));
                return;
            }

        /* Output first line of description */
        /* Find good place for next line break */
        for (j = l - n - 3; j > 0; j--)
        {
            if ((description[j] == '-') && (description[j - 1] != ' '))
            {
                /* Break after hyphen */
                j++;
                break;
            }
            if ((description[j] == ' ') && (description[j + 1] != '-'))
            {
                /* Break at space, unless it borders a dash */
                while ((j > 0) && (description[j] == ' '))
                    j--;
                if ((j > 0) && (description[j - 1] != '-'))
                    break;
            }
        }
        if (j == 0)
        {
            /* No good line break found - look for acceptable line break */
            for (j = l - n - 4; (j > 0) && (description[j] != ' '); j--)
                ;
            while (description[j] == ' ')
                j--;
            if (j == 0)
                j = l - n - 4;
        }
        fwrite(description, sizeof(char), j + 1, fp);
        fputc('\n', fp);
        description += j + 1;
        /* Output remainder of description */
        pretty_print(fp, description, l, (i < 0 ? n + 2 : i));
    }
}

static void _print_usage(FILE *f, char *name)
{
    fprintf(f, "Usage: %s [options] < [input]\n", name);
    pretty_print(f, "The program reads data from the input file specified and constructs history by threading a sequence at a time.", 70, 0);
    fprintf(f, "Legal options are:\n");
    print_option(f, "-M[x]", "Specify cost of a recurrent mutation (default: x = 0.9).", 70, -1);
    print_option(f, "-B[x]", "Specify cost of a back mutation (default: x = 1.0).", 70, -1);
    print_option(f, "-R[x]", "Specify cost of a single recombination (default: x = 1.1).", 70, -1);
    print_option(f, "-C[x]", "Specify cost of a multiple cross-over recombination (default: x = 1.1).", 70, -1);
    print_option(f, "-V[x]", "level of verbosity", 70, -1);
    print_option(f, "-d[name]", "Output ancestral recombination graph of minimum recombination history in dot format to file name.", 70, -1);
    print_option(f, "-g[name]", "Output ancestral recombination graph of minimum recombination history in GDL format to file name.", 70, -1);
    print_option(f, "-j[name]", "Output ancestral recombination graph of minimum recombination history in GML format to file name.", 70, -1);
    print_option(f, "-e[x]", "How to label edges in output. 0: not at all, 1: recurrent and back mutations, 2: everything", 70, -1);
    print_option(f, "-s", "Label nodes in ancestral recombination graphs with mutations of that node.", 70, -1);
    print_option(f, "-l", "Input data has labels for the samples.", 70, -1);
    print_option(f, "-i", "Label sequences by index, e.g. the first sequence in the data file would be assigned '1' as id.", 70, -1);
    print_option(f, "-r[x]", "Input data has x root(s) given. Otherwise root is all zero. More than one root has odd behaviour TODO:", 70, -1);
    print_option(f, "-f[x]", "Find root for given sequences, using strategy x", 70, -1);
    print_option(f, "-Q[x]", "Sets the number of runs.", 70, -1);
    print_option(f, "-S[x]", "Sets the random seed.", 70, -1);
    print_option(f, "-X[x]", "Provide an upper bound x on the number of recombinations needed for the input dataset.", 70, -1);
    print_option(f, "-Y[x]", "Provide an upper bound x on the number of recurrent mutations needed for the input dataset.", 70, -1);
    print_option(f, "-Z[x]", "Provide an upper bound x on the number of back mutations needed for the input dataset.", 70, -1);
    print_option(f, "-T[x]", "Run type (used on multiruns). 0: simple random, 1: deal with problem sites first, 2: deal with problem sites last", 70, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 70, -1);
}

Genes read_input(std::istream &in, bool use_labels, bool generate_ids)
{
    Genes genes;
    std::string line;
    Gene g;
    int seq_index = 0;

    if (!use_labels)
    {
        int i = 1;
        while (std::getline(in, line))
        {
            if (generate_ids)
                g.label = std::to_string(i);
            else
                g.label = "";
            g.mutations.clear();
            seq_index = 0;
            for (auto &c : line)
            {
                if (c == '1')
                {
                    g.mutations.push_back(seq_index);
                }
                seq_index += 1;
            }
            genes.genes.push_back(g); // Makes copy
            i += 1;
        }
    }
    else
    {

        bool is_first_line = true;
        while (std::getline(in, line))
        {
            if (line[0] == '>')
            {
                if (!is_first_line) //
                {
                    genes.genes.push_back(g); // This makes a copy
                    g.label = "";
                    g.mutations.clear();
                    seq_index = 0;
                }
                g.label = line.substr(1);
            }
            else
            {
                for (auto &c : line)
                {
                    if (c == '1')
                    {
                        g.mutations.push_back(seq_index);
                    }
                    seq_index += 1;
                }
            }
            is_first_line = false;
        }
        genes.genes.push_back(g); // Make sure final sequence added
    }

    return std::move(genes);
}

/* flips all sites which are mutated in root */
void set_genes_relative_to_root(Genes &genes, const Gene &root)
{
    for (Gene &g : genes.genes)
    {
        g.mutations = vector_symmetric_difference(g.mutations, root.mutations);
    }
}

void set_arg_relative_to_root(ARG &arg, const Gene &root)
{
    for (auto &node : arg.nodes)
    {
        node->mutations = vector_symmetric_difference(node->mutations, root.mutations);
    }
}

int main(int argc, char **argv)
{
    FILE *fp;

    int how_to_label_edges = 0; // 0 is not at all, 1 is bm and rms, 2 is everything.
    bool do_label_node_mutations = false;
    bool do_generate_ids = false;

    bool root_given = false;
    int number_roots_given = 0;
    int find_root_strategy = -1; // -1 means root is given/assumed to be all zero

    int run_seed = 0;
    int num_runs = 0;
    int how_verbose = 0;
    bool label_sequences = false;

    float cost_rm = 0.9;
    float cost_bm = 1.0;
    float cost_recomb = 1.1;
    float cost_multirecomb = 1.1;

    // Negative values mean unlimited
    int recomb_max = -1, rm_max = -1, bm_max = -1;

    int run_type = 0;

    std::vector<std::string> dot_files = {};
    std::vector<std::string> gml_files = {};
    std::vector<std::string> gdl_files = {};

#define KWARG_OPTIONS "M:B:R:C:V:d::g::j::e:slir:f:Q:S:X:Y:Z:T:hH?"

    /* Parse command line options */
    int i;
    while ((i = getopt(argc, argv, KWARG_OPTIONS)) >= 0)
    {
        switch (i)
        {
        case 'M':
            cost_rm = std::stof(optarg);
            if (errno != 0)
            {
                std::cerr << "Some error occured with cost_rm input.\n";
                exit(1);
            }
            if (cost_rm < 0 && cost_rm != 0)
            {
                std::cerr << "negative value (normally -1) means recurrent mutations not allowed.\n";
            }
            break;
        case 'B':
            cost_bm = std::stof(optarg);
            if (errno != 0)
            {
                std::cerr << "Some error occured with cost_bm input.\n";
                exit(1);
            }
            if (cost_bm < 0 && cost_bm != 0)
            {
                std::cerr << "negative value (normally -1) means back mutations not allowed.\n";
            }
            break;
        case 'R':
            cost_recomb = std::stof(optarg);
            if (errno != 0)
            {
                std::cerr << "Some error occured with cost_recomb input.\n";
                exit(1);
            }
            if (cost_recomb < 0 && cost_recomb != 0)
            {
                std::cerr << "negative value (normally -1) means recombinations not allowed.\n";
            }
            break;
        case 'C':
            cost_multirecomb = std::stof(optarg);
            if (errno != 0)
            {
                std::cerr << "Some error occured with cost_multirecomb input.\n";
                exit(1);
            }
            if (cost_multirecomb < 0 && cost_multirecomb != 0)
            {
                std::cerr << "negative value (normally -1) means multiple recombinations not allowed.\n";
            }
            break;
        case 'V':
            how_verbose = std::stoi(optarg);
            if (errno != 0 || (how_verbose > 3 && how_verbose < 0))
            {
                std::cerr << "Verbosity input should be between 0 and 3 inclusive.\n";
                exit(1);
            }
            break;
        case 'd':
            /* Output ancestral recombination graph of history leading to
             * minimum number of recombinations in dot format.
             */
            /* Was a file name specified? */
            if (optarg != 0)
            {
                if (optarg[0] == '-')
                {
                    std::cerr << "Option -d requires an output file.\n";
                    exit(1);
                }
                /* Check whether file can be written before initiating computation */
                if ((fp = fopen(optarg, "w")) == NULL)
                {
                    std::cerr << "Could not open file " << optarg << "for output.\n";
                    exit(1);
                }
                fclose(fp);
                dot_files.push_back(optarg);
            }
            else
                std::cerr << "Please specify file\n";
            break;
        case 'g':
            /* Output ancestral recombination graph of history leading to
             * minimum number of recombinations in gdl format.
             */
            /* Was a file name specified? */
            if (optarg != 0)
            {
                if (optarg[0] == '-')
                {
                    std::cerr << "Option -g requires an output file.\n";
                    exit(1);
                }
                /* Check whether file can be written before initiating computation */
                if ((fp = fopen(optarg, "w")) == NULL)
                {
                    std::cerr << "Could not open file " << optarg << "for output.\n";
                    exit(1);
                }
                fclose(fp);
                gdl_files.push_back(optarg);
            }
            else
                std::cerr << "Please specify file\n";
            break;
        case 'j':
            /* Output ancestral recombination graph of history leading to
             * minimum number of recombinations in gml format.
             */
            /* Was a file name specified? */
            if (optarg != 0)
            {
                if (optarg[0] == '-')
                {
                    std::cerr << "Option -j requires an output file.\n";
                    exit(1);
                }
                /* Check whether file can be written before initiating computation */
                if ((fp = fopen(optarg, "w")) == NULL)
                {
                    std::cerr << "Could not open file " << optarg << "for output.\n";
                    exit(1);
                }
                fclose(fp);
                gml_files.push_back(optarg);
            }
            else
                std::cerr << "Please specify file\n";
            break;
        case 'e':
            how_to_label_edges = std::stoi(optarg);
            if (errno != 0 || how_to_label_edges < 0 || how_to_label_edges > 2)
            {
                std::cerr << "how_to_label_edges should be between 0 and 2.\n";
                exit(1);
            }
            break;
        case 's':
            do_label_node_mutations = true;
            break;
        case 'l':
            label_sequences = true;
            break;
        case 'i':
            do_generate_ids = true;
            break;
        case 'r':
            number_roots_given = std::stoi(optarg);
            if (errno != 0 || number_roots_given < 0)
            {
                std::cerr << "Number of roots given should be a non-negative integer.\n";
                exit(1);
            }
            if (number_roots_given > 0)
                root_given = true;
            break;
        case 'f':
            find_root_strategy = std::stoi(optarg);
            if (errno != 0 || find_root_strategy < 0)
            {
                std::cerr << "Number of root strategy should be a non-negative integer.\n";
                exit(1);
            }
            break;
        case 'Q':
            num_runs = std::stoi(optarg);
            if (errno != 0 || num_runs < 0)
            {
                std::cerr << "Number of iterations should be a positive integer.\n";
                exit(1);
            }
            break;
        case 'S':
            run_seed = std::stoi(optarg);
            if (errno != 0 || run_seed < 0)
            {
                std::cerr << "Seed input should be a positive integer.\n";
                exit(1);
            }
            break;
        case 'X':
            recomb_max = std::stoi(optarg);
            if (errno != 0 || recomb_max < 0)
            {
                std::cerr << "Upper bound on number of recombinations should be a positive integer.\n";
                exit(1);
            }
            break;
        case 'Y':
            rm_max = std::stoi(optarg);
            if (errno != 0 || rm_max < 0)
            {
                std::cerr << "Upper bound on number of recurrent mutations should be a positive integer.\n";
                exit(1);
            }
            break;
        case 'Z':
            bm_max = std::stoi(optarg);
            if (errno != 0 || bm_max < 0)
            {
                std::cerr << "Upper bound on number of recurrent mutations should be a positive integer.\n";
                exit(1);
            }
            break;
        case 'T':
            run_type = std::stoi(optarg);
            if (errno != 0 || run_type < 0)
            {
                std::cerr << "Run type should be a non-negative integer.\n";
                exit(1);
            }
            break;
        case 'h':
        case 'H':
        case '?':
            _print_usage(stdout, argv[0]);
            exit(0);
        case ':':
            _print_usage(stderr, argv[0]);
            exit(1);
        }
    }

    if (number_roots_given > 0 && find_root_strategy >= 0)
    {
        std::cerr << "Can't have roots given whilst also wanting to find roots.\n";
        exit(1);
    }
    if (do_generate_ids && label_sequences)
    {
        std::cerr << "Can't generate labels if labels are also meant to be read.\n";
        exit(1);
    }

    if (how_verbose >= 1)
        std::cout << "parsed inputs\n";

    // TODO: Read input from file or stdin
    Genes genes = read_input(std::cin, label_sequences, do_generate_ids);

    if (how_verbose >= 1)
        std::cout << "read input\n";

    Gene root;
    if (root_given && number_roots_given == 1)
    {
        root = genes.genes[0];
        set_genes_relative_to_root(genes, root);
    }
    else
        root.label = "root";

    ARG arg;
    if (num_runs == 0)
    {
        arg = build_arg(genes, root_given, how_verbose, cost_rm, cost_bm, cost_recomb, recomb_max, rm_max, bm_max);
    }
    else
    {
        switch (run_type)
        {
        case 0:
            arg = build_arg_multi_random_runs(num_runs, run_seed, genes, root_given, how_verbose, cost_rm, cost_bm, cost_recomb, recomb_max, rm_max, bm_max);
            break;
        case 1:
            arg = build_arg_multi_smart_runs(num_runs, run_seed, true, genes, root_given, how_verbose, cost_rm, cost_bm, cost_recomb, recomb_max, rm_max, bm_max);
            break;
        case 2:
            arg = build_arg_multi_smart_runs(num_runs, run_seed, false, genes, root_given, how_verbose, cost_rm, cost_bm, cost_recomb, recomb_max, rm_max, bm_max);
            break;
        }
    }


    if (root_given && number_roots_given == 1)
    {
        set_arg_relative_to_root(arg, root);
    }

    /* Output ARG in dot format */
    for (auto dot_file : dot_files)
    {
        fp = fopen(dot_file.c_str(), "w");
        auto label_format = LABEL_LABEL;
        if (do_label_node_mutations)
            label_format = LABEL_BOTH;
        arg_output(arg, genes, fp, ARGDOT, how_to_label_edges, label_format);
        fclose(fp);
    }

    /* Output GML in dot format */
    for (auto gml_file : gml_files)
    {
        fp = fopen(gml_file.c_str(), "w");
        auto label_format = LABEL_LABEL;
        if (do_label_node_mutations)
            label_format = LABEL_BOTH;
        arg_output(arg, genes, fp, ARGGML, how_to_label_edges, label_format);
        fclose(fp);
    }

    /* Output GDL in dot format */
    for (auto gdl_file : gdl_files)
    {
        fp = fopen(gdl_file.c_str(), "w");
        auto label_format = LABEL_LABEL;
        if (do_label_node_mutations)
            label_format = LABEL_BOTH;
        arg_output(arg, genes, fp, ARGGDL, how_to_label_edges, label_format);
        fclose(fp);
    }

    return 0;
}
