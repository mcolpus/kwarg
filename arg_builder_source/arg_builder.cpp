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
    print_option(f, "-L[x]", "Provide an optional label x (should be an integer) to print at the start of each line.", 70, -1);
    print_option(f, "-S[x]", "Specify cost of a sequencing error (default: x = 0.5).", 70, -1);
    print_option(f, "-M[x]", "Specify cost of a recurrent mutation (default: x = 0.9).", 70, -1);
    print_option(f, "-R[x]", "Specify cost of a single recombination (default: x = 1.0).", 70, -1);
    print_option(f, "-C[x]", "Specify cost of two recombinations immediately following each other (default: x = 2.0).", 70, -1);
    print_option(f, "-V[x]", "level of verbosity", 70, -1);
    print_option(f, "-d[name]", "Output ancestral recombination graph of minimum recombination history in dot format to file name.", 70, -1);
    print_option(f, "-g[name]", "Output ancestral recombination graph of minimum recombination history in GDL format to file name.", 70, -1);
    print_option(f, "-j[name]", "Output ancestral recombination graph of minimum recombination history in GML format to file name.", 70, -1);
    print_option(f, "-i", "Sequences not having a sequence id in the data file are assigned their index in the data file as id, e.g. the first sequence in the data file would be assigned '1' as id.", 70, -1);
    print_option(f, "-e", "Label edges in ancestral recombination graphs with the sites undergoing mutation along the edge.", 70, -1);
    print_option(f, "-o", "Assume input data is in own format. Default is to first try to parse data in own format, and if that fails to try to parse it in fasta format. Specifying this option, no attempt will be made to try to parse the data in fasta format.", 70, -1);
    print_option(f, "-f", "Assume input data is in fasta format. No attempt will be made to try to parse the data in own format. Note that the -o and the -f options override each other, so only the last one occurring in the command line will have an effect.", 70, -1);
    print_option(f, "-a", "Assume input consists of amino acid (protein) sequences using the one letter amino acid alphabet; anything not in the amino acid one letter alphabet is treated as an unresolved site (default is to assume sequences in binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known (see option -k), the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 70, -1);
    print_option(f, "-n", "Assume input consists of nucleic sequences; anything but a/c/g/t/u is considered an unresolved site (default is to assume binary, i.e. 0/1, format where anything but a 0 or a 1 is considered an unresolved site). If the most recent common ancestor is assumed known (see option -k), the first sequence in the input data is considered to specify the wild type of each site and is not included as part of the sample set.", 70, -1);
    print_option(f, "-l", "Input data has labels for the samples.", 70, -1);
    print_option(f, "-Q[x]", "Sets the number of runs.", 70, -1);
    print_option(f, "-Z[x]", "Sets the random seed z (only one run is made in this case).", 70, -1);
    print_option(f, "-X[x]", "Provide an upper bound x on the number of recombinations needed for the input dataset (solutions with more than x recombinations will be abandoned).", 70, -1);
    print_option(f, "-Y[x]", "Provide an upper bound x on the number of recurrent needed for the input dataset (solutions with more than x recurrent mutations will be abandoned).", 70, -1);
    print_option(f, "-h, -H -?", "Print this information and stop.", 70, -1);
}

Genes read_input(std::istream &in, bool use_labels)
{
    Genes genes;

    Gene g;
    int seq_index = 0;

    std::string line;

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
    FILE *print_progress = stdout;
    FILE *fp;

    Gene_Format format = GENE_ANY;
    Gene_SeqType seqtype = GENE_BINARY;

    int i;
    char *token;
    int run_reference = -1;
    bool do_label_edges = false;
    bool do_generate_ids = false;

    int run_seed = 0;
    int num_runs = 0;
    int how_verbose = 0;
    bool label_sequences = false;

    int rec_max = 0, rm_max = 0;

    std::vector<std::string> dot_files = {};
    std::vector<std::string> gml_files = {};
    std::vector<std::string> gdl_files = {};

#define KWARG_OPTIONS "L:S:M:R:C:V:d::g::j::ieofanlQ:Z:X:Y:hH?"

    /* Parse command line options */
    while ((i = getopt(argc, argv, KWARG_OPTIONS)) >= 0)
    {
        switch (i)
        {
        case 'L':
            run_reference = std::stoi(optarg);
            if (errno != 0 || run_reference < 0)
            {
                fprintf(stderr, "Reference should be a positive integer.\n");
                exit(1);
            }
            break;
        case 'S':
            if (optarg != 0)
            {
                fprintf(stderr, "Not yet implemented.\n");
            }
            break;
        case 'M':
            if (optarg != 0)
            {
                fprintf(stderr, "Not yet implemented.\n");
            }
            break;
        case 'R':
            if (optarg != 0)
            {
                fprintf(stderr, "Not yet implemented.\n");
            }
            break;
        case 'C':
            if (optarg != 0)
            {
                fprintf(stderr, "Not yet implemented.\n");
            }
            break;
        case 'V':
            how_verbose = std::stoi(optarg);
            if (errno != 0 || (how_verbose > 2 && how_verbose < 0))
            {
                fprintf(stderr, "Verbosity input should be 0, 1 or 2.\n");
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
                    fprintf(stderr, "Option -d requires an output file.\n");
                    exit(1);
                }
                /* Check whether file can be written before initiating computation */
                if ((fp = fopen(optarg, "w")) == NULL)
                {
                    fprintf(stderr, "Could not open file %s for output\n", optarg);
                    exit(1);
                }
                fclose(fp);
                dot_files.push_back(optarg);
            }
            else
                fprintf(stderr, "Please specify file\n"); // TODO: print to stdout
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
                    fprintf(stderr, "Option -g requires an output file.\n");
                    exit(1);
                }
                /* Check whether file can be written before initiating computation */
                if ((fp = fopen(optarg, "w")) == NULL)
                {
                    fprintf(stderr, "Could not open file %s for output\n", optarg);
                    exit(1);
                }
                fclose(fp);
                gdl_files.push_back(optarg);
            }
            else
                fprintf(stderr, "Please specify file\n"); // TODO: print to stdout
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
                    fprintf(stderr, "Option -j requires an output file.\n");
                    exit(1);
                }
                /* Check whether file can be written before initiating computation */
                if ((fp = fopen(optarg, "w")) == NULL)
                {
                    fprintf(stderr, "Could not open file %s for output\n", optarg);
                    exit(1);
                }
                fclose(fp);
                gml_files.push_back(optarg);
            }
            else
                fprintf(stderr, "Please specify file\n"); // TODO: print to stdout
            break;
        case 'i':
            do_generate_ids = true;
            break;
        case 'e':
            do_label_edges = true;
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
        case 'Q':
            num_runs = std::stoi(optarg);
            if (errno != 0 || num_runs < 0)
            {
                fprintf(stderr, "Number of iterations should be a positive integer.\n");
                exit(1);
            }
            break;
        case 'Z':
            run_seed = std::stoi(optarg);
            if (errno != 0 || run_seed <= 0)
            {
                fprintf(stderr, "Seed input should be a positive integer.\n");
                exit(1);
            }
            break;
        case 'X':
            rec_max = std::stoi(optarg);
            if (errno != 0 || rec_max < 0)
            {
                fprintf(stderr, "Upper bound on number of recombinations should be a positive integer.\n");
                exit(1);
            }
            break;
        case 'Y':
            rm_max = std::stoi(optarg);
            if (errno != 0 || rm_max < 0)
            {
                fprintf(stderr, "Upper bound on number of recurrent mutations should be a positive integer.\n");
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

    fprintf(print_progress, "parsed inputs\n");

    // TODO: Read input from file or stdin
    Genes genes = read_input(std::cin, label_sequences);

    fprintf(print_progress, "read input\n");

    build_arg(genes, print_progress);

    return 0;
}
