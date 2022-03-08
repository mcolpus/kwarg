from cmath import inf
import re
import numpy as np
import sys
import getopt
import matplotlib as mpl
from matplotlib import pyplot as plt
from typing import List, Dict, Tuple


np.set_printoptions(edgeitems=30, linewidth=100000,
                    formatter=dict(float=lambda x: "%.3g" % x))

plt.rc('font', size=12) 

# Reads from csv containing many runs and create a grid for optimal args


def read_args(inputfile):
    arg_runs = []

    file = open(inputfile, "r")
    for line in file.readlines():
        if len(line) == 0:
            continue
        elif line[0] == 'r':
            continue  # this is the header which starts recombination

        values = line.split(',')
        recombinations = int(values[0])
        back_mutations = int(values[1]) # On kwarg this is sequential errors
        recurrent_mutations = int(values[2])
        if (recombinations != -1):
            arg_runs.append((recombinations, back_mutations, recurrent_mutations))

    file.close()

    return arg_runs


def runs_to_dict(arg_runs):
    fewest_rms = {}

    max_recombs = 0
    max_bms = 0

    for (recombs, bms, rms) in arg_runs:
        max_recombs = max(max_recombs, recombs)
        max_bms = max(max_bms, bms)

        if (recombs, bms) in fewest_rms:
            current_rms = fewest_rms[(recombs, bms)]
            if (rms < current_rms):
                fewest_rms[(recombs, bms)] = rms
        else:
            fewest_rms[(recombs, bms)] = rms

    return (fewest_rms, max_recombs, max_bms)


def calculate_optimal_grid(fewest_rms, max_recombs, max_bms):
    grid = np.full((max_recombs+1, max_bms+1), sys.maxsize, dtype=int)

    # Fill out grid with known args first
    for ((recombs, bms), rms) in fewest_rms.items():
        grid[recombs, bms] = rms

    # Now fill in grid
    for i in range(max_recombs+1):
        for j in range(max_bms+1):
            if i > 0 and j > 0:
                grid[i, j] = min(grid[i, j], grid[i-1, j], grid[i, j-1])
            elif i > 0:
                grid[i, j] = min(grid[i, j], grid[i-1, j])
            elif j > 0:
                grid[i, j] = min(grid[i, j], grid[i, j-1])
    
    for i in range(max_recombs+1):
        for j in range(max_bms+1):
            if grid[i, j] == sys.maxsize:
                grid[i, j] = -1

    return grid


def graph_grid(grid, outputfile):
    fig = plt.figure(1)

    ax1 = fig.add_axes([0.15, 0.12, 0.65, 0.75])

    n = np.max(grid)
    colours = ['green', 'blue', 'red']
    nodes = [0.0, 0.5, 1.0]
    bounds = range(0,n+2)

    if (np.min(grid) == -1):
        colours = ['black', 'green', 'blue', 'red']
        nodes = [0.0, 1.0/float(n), 0.5, 1.0]
        bounds = range(-1,n+2)

    cmap = mpl.colors.LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colours)), N = len(bounds) - 1)

    img = plt.imshow(grid.transpose(), interpolation='nearest',
                     cmap=cmap, aspect='auto', origin='lower')

    plt.title(
        'Fewest recurrent-mutations required \n given recombination and back-mutations')
    plt.xlabel('Recombinations')
    plt.ylabel('Back-mutations')

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax2 = fig.add_axes([0.85, 0.12, 0.05, 0.75])
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=ax2, orientation='vertical',
                label="recurrent mutations")

    fig.savefig(outputfile + "_grid.png")


def graph_all_recombs_vs_rare_mutations(arg_runs, outputfile):
    fig = plt.figure(2)

    recombinations = []
    rare_muts = []
    max_recombs = 0
    for (recombs, bms, rms) in arg_runs:
        recombinations.append(recombs)
        rare_muts.append(bms + rms)
        max_recombs = max(max_recombs, recombs)

    fewest_rare_muts = [sys.maxsize for i in range(max_recombs+1)]
    for (recombs, bms, rms) in arg_runs:
        fewest_rare_muts[recombs] = min(fewest_rare_muts[recombs], bms+rms)

    for i in range(1, max_recombs+1):
        fewest_rare_muts[i] = min(fewest_rare_muts[i-1], fewest_rare_muts[i])

    plt.plot(range(max_recombs+1), fewest_rare_muts,
             color='red', label="optimal")
    plt.scatter(recombinations, rare_muts, color='blue',
                label='all networks created')

    plt.title('Recombinations vs Recurrent mutations \n for all args created')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent mutations (all kinds)')
    plt.legend(loc="upper right")

    # plt.show()

    fig.savefig(outputfile + "_rare_muts.png")

def graph_comparison(arg_runs_list, labels, outputfile):
    fig = plt.figure(3)

    for arg_runs, label in zip(arg_runs_list, labels):
        recombinations = []
        rare_muts = []
        max_recombs = 0
        for (recombs, bms, rms) in arg_runs:
            recombinations.append(recombs)
            rare_muts.append(bms + rms)
            max_recombs = max(max_recombs, recombs)
        
        fewest_rare_muts = [sys.maxsize for i in range(max_recombs+1)]
        for (recombs, bms, rms) in arg_runs:
            fewest_rare_muts[recombs] = min(fewest_rare_muts[recombs], bms+rms)
        for i in range(1, max_recombs+1):
            fewest_rare_muts[i] = min(fewest_rare_muts[i-1], fewest_rare_muts[i])
        
        plt.plot(range(max_recombs+1), fewest_rare_muts, label=label, alpha=0.7)

    plt.title('Recombinations vs Recurrent Mutations')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent mutations (all kinds)')
    plt.legend(loc="upper right")

    fig.savefig(outputfile + ".png")


def main(argv):
    inputfile = ''
    outputfile = ''
    labels_text = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:l:", ["ifile=", "ofile=", "labels="])
    except getopt.GetoptError:
        print('calculate_optimal_arg_grid.py -i <input file> -o <output file> -l <labels>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('calculate_optimal_arg_grid.py -i <input file> -o <output file> -l <labels>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-l", "--labels"):
            labels_text = arg

    print('Input file is: ', inputfile)
    print('Output file is: ', outputfile)
    print('Labels are: ', labels_text)

    if ',' in inputfile:
        # Perform a comparison
        inputfiles = inputfile.split(',')
        arg_runs_list = [read_args(file) for file in inputfiles]
        graph_comparison(arg_runs_list, labels_text.split(','), outputfile)
    else:
        arg_runs = read_args(inputfile)
        (fewest_rms, max_recombs, max_bms) = runs_to_dict(arg_runs)
        grid = calculate_optimal_grid(fewest_rms, max_recombs, max_bms)

        graph_grid(grid, outputfile)

        graph_all_recombs_vs_rare_mutations(arg_runs, outputfile)

        print(max_recombs)
        print(max_bms)
        print(fewest_rms)
        print(grid)


if __name__ == "__main__":
    main(sys.argv[1:])
