from cmath import inf
import re
import numpy as np
import sys
import getopt
import matplotlib as mpl
from matplotlib import pyplot
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
        back_mutations = int(values[1])
        recurrent_mutations = int(values[2])
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

    return grid


def graph_grid(grid, outputfile):
    fig = pyplot.figure(1)

    cmap2 = mpl.colors.LinearSegmentedColormap.from_list('my_colormap',
                                                         ['green', 'blue', 'red'],
                                                         256)

    img2 = pyplot.imshow(grid, interpolation='nearest',
                         cmap=cmap2, aspect='auto',
                         origin='lower')

    pyplot.colorbar(img2, cmap=cmap2)

    pyplot.show()

    fig.savefig(outputfile)


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

    arg_runs = read_args(inputfile)
    (fewest_rms, max_recombs, max_bms) = runs_to_dict(arg_runs)
    grid = calculate_optimal_grid(fewest_rms, max_recombs, max_bms)

    graph_grid(grid, outputfile)

    print(max_recombs)
    print(max_bms)
    print(fewest_rms)
    print(grid)


if __name__ == "__main__":
    main(sys.argv[1:])
