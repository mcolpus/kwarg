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

    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                        ['green', 'blue', 'red'],
                                                        np.max(grid) - np.min(grid) + 1)

    img = plt.imshow(grid.transpose(), interpolation='nearest',
                     cmap=cmap, aspect='auto', origin='lower')

    plt.colorbar(img, cmap=cmap)

    plt.title(
        'Fewest recurrent-mutations required \n given recombination and back-mutations')
    plt.xlabel('Recombinations')
    plt.ylabel('Back-mutations')

    plt.show()

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

    plt.title('Recombinations vs Rare Mutations \n for all args created')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent + back mutations')
    plt.legend(loc="upper right")

    # plt.show()

    fig.savefig(outputfile + "_rare_muts.png")

def graph_comparison(arg_runs1, arg_runs2, outputfile):
    fig = plt.figure(3)

    recombinations1 = []
    rare_muts1 = []
    max_recombs = 0
    for (recombs, bms, rms) in arg_runs1:
        recombinations1.append(recombs)
        rare_muts1.append(bms + rms)
        max_recombs = max(max_recombs, recombs)

    recombinations2 = []
    rare_muts2 = []
    for (recombs, bms, rms) in arg_runs2:
        recombinations2.append(recombs)
        rare_muts2.append(bms + rms)
        max_recombs = max(max_recombs, recombs)

    fewest_rare_muts1 = [sys.maxsize for i in range(max_recombs+1)]
    for (recombs, bms, rms) in arg_runs1:
        fewest_rare_muts1[recombs] = min(fewest_rare_muts1[recombs], bms+rms)
    for i in range(1, max_recombs+1):
        fewest_rare_muts1[i] = min(fewest_rare_muts1[i-1], fewest_rare_muts1[i])
    
    fewest_rare_muts2 = [sys.maxsize for i in range(max_recombs+1)]
    for (recombs, bms, rms) in arg_runs2:
        fewest_rare_muts2[recombs] = min(fewest_rare_muts2[recombs], bms+rms)
    for i in range(1, max_recombs+1):
        fewest_rare_muts2[i] = min(fewest_rare_muts2[i-1], fewest_rare_muts2[i])

    plt.plot(range(max_recombs+1), fewest_rare_muts1,
             color='red', label="optimal 1")
    plt.plot(range(max_recombs+1), fewest_rare_muts2,
             color='blue', label="optimal 2")
    plt.scatter(recombinations1, rare_muts1, color='green', alpha=0.5,
                label='all networks created 1')
    plt.scatter(recombinations2, rare_muts2, color='yellow', alpha=0.5,
                label='all networks created 2')

    plt.title('Recombinations vs Rare Mutations \n comparison')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent + back mutations')
    plt.legend(loc="upper right")

    # plt.show()

    fig.savefig(outputfile + ".png")


def main(argv):
    inputfile = ''
    outputfile = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('clean_binary_sequences.py -i <input file> -o <output file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('clean_binary_sequences.py -i <input file> -o <output file>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print('Input file is: ', inputfile)
    print('Output file is: ', outputfile)

    if ',' in inputfile:
        inputfile1, inputfile2 = inputfile.split(',')
        arg_runs1 = read_args(inputfile1)
        arg_runs2 = read_args(inputfile2)
        graph_comparison(arg_runs1, arg_runs2, outputfile)
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
