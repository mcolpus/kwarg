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
    print_option(f, "-V[x]", "level of verbosity", 70, -1);
    print_option(f, "-d[name]", "Output ancestral recombination graph of minimum recombination history in dot format to file name.", 70, -1);
    print_option(f, "-g[name]", "Output ancestral recombination graph of minimum recombination history in GDL format to file name.", 70, -1);
    print_option(f, "-j[name]", "Output ancestral recombination graph of minimum recombination history in GML format to file name.", 70, -1);
    print_option(f, "-i", "Sequences not having a sequence id in the data file are assigned their index in the data file as id, e.g. the first sequence in the data file would be assigned '1' as id.", 70, -1);
    print_option(f, "-e", "Label edges in ancestral recombination graphs with the sites undergoing mutation along the edge.", 70, -1);
    print_option(f, "-s", "Label nodes in ancestral recombination graphs with mutations of that node.", 70, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 70, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 70, -1);
    print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known (see option -k), the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 70, -1);
    print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known (see option -k), the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 70, -1);
    print_option(f, "-l", "Input data has labels for the samples.", 70, -1);
    print_option(f, "-r", "Input data has root as first sample. If not than root is assumed to be all zero", 70, -1);
    print_option(f, "-Q[x]", "Sets the number of runs.", 70, -1);
    print_option(f, "-S[x]", "Sets the random seed.", 70, -1);
    print_option(f, "-X[x]", "Provide an upper bound x on the number of recombinations needed for the input dataset.", 70, -1);
    print_option(f, "-Y[x]", "Provide an upper bound x on the number of recurrent mutations needed for the input dataset.", 70, -1);
    print_option(f, "-Z[x]", "Provide an upper bound x on the number of back mutations needed for the input dataset.", 70, -1);
    print_option(f, "-T[x]", "Run type (used on multiruns). 0: simple random, 1: deal with problem sites first, 2: deal with problem sites last", 70, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 70, -1);
}

Genes read_input(std::istream &in, bool use_labels)
{
    Genes genes;

    // std::string root_line;
    // // If root given than it is first line
    // if (root_given)
    // {
    //     if (use_labels)
    //     {
    //         std::getline(in, root_line);
    //         if (root_line[0] != '>')
    //             std::cerr << "root doesn't have label";
    //     }
    // }

    std::string line;
    Gene g;
    int seq_index = 0;


    if (!use_labels)
    {
        while (std::getline(in, line))
        {
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
    std::cout << "done\n";

    return std::move(genes);
}

int main(int argc, char **argv)
{
    FILE *fp;

    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;

    int i;
    char *token;
    int run_reference = -1;
    bool do_label_edges = false;
    bool do_label_node_mutations = false;
    bool do_generate_ids = false;
    bool root_given = false;

    int run_seed = 0;
    int num_runs = 0;
    int how_verbose = 0;
    bool label_sequences = false;

    float cost_rm = 0.9;
    float cost_bm = 1.0;
    float cost_recomb = 1.1;
    
    // Negative values mean unlimited
    int recomb_max = -1, rm_max = -1, bm_max = -1;

    int run_type = 0;

    std::vector<std::string> dot_files = {};
    std::vector<std::string> gml_files = {};
    std::vector<std::string> gdl_files = {};

#define KWARG_OPTIONS "L:M:B:R:V:d::g::j::iesofanlrQ:S:X:Y:Z:T:hH?"

    /* Parse command line options */
    while ((i = getopt(argc, argv, KWARG_OPTIONS)) >= 0)
    {
        switch (i)
        {
        case 'L':
            run_reference = std::stoi(optarg);
            if (errno != 0 || run_reference < 0)
            {
                std::cerr << "Reference should be a positive integer.\n";
                exit(1);
            }
            break;
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
        case 'i':
            do_generate_ids = true;
            break;
        case 'e':
            do_label_edges = true;
            break;
        case 's':
            do_label_node_mutations = true;
            break;
        case 'o':
            format = GENE_BEAGLE;
            break;
        case 'f':
            format = GENE_FASTA;
            break;
        case 'a':
            seqtype = GENE_AMINO;
            break;
        case 'n':
            seqtype = GENE_NUCLEIC;
            break;
        case 'l':
            label_sequences = true;
            break;
        case 'r':
            root_given = true;
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

    if (how_verbose >= 1)
        std::cout << "parsed inputs\n";

    // TODO: Read input from file or stdin
    Genes genes = read_input(std::cin, label_sequences);

    if (how_verbose >= 1)
        std::cout << "read input\n";

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

    /* Output ARG in dot format */
    for (auto dot_file : dot_files)
    {
        fp = fopen(dot_file.c_str(), "w");
        auto label_format = LABEL_LABEL;
        if (do_label_node_mutations)
            label_format = LABEL_BOTH;
        arg_output(arg, genes, fp, ARGDOT, do_label_edges, label_format);
        fclose(fp);
    }

    /* Output GML in dot format */
    for (auto gml_file : gml_files)
    {
        fp = fopen(gml_file.c_str(), "w");
        auto label_format = LABEL_LABEL;
        if (do_label_node_mutations)
            label_format = LABEL_BOTH;
        arg_output(arg, genes, fp, ARGGML, do_label_edges, label_format);
        fclose(fp);
    }

    /* Output GDL in dot format */
    for (auto gdl_file : gdl_files)
    {
        fp = fopen(gdl_file.c_str(), "w");
        auto label_format = LABEL_LABEL;
        if (do_label_node_mutations)
            label_format = LABEL_BOTH;
        arg_output(arg, genes, fp, ARGGDL, do_label_edges, label_format);
        fclose(fp);
    }

    return 0;
}
