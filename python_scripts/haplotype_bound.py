from re import sub
from xmlrpc.client import Unmarshaller
import numpy as np
import sys
import getopt
import itertools
from itertools import chain , combinations

genes = np.zeros(0)
num_sequences = 0
num_cols = 0


def read_sequences(inputfile):
    global genes
    global num_sequences
    global num_cols

    sequences = []
    sequence_labels = []

    file = open(inputfile, "r")
    seq = ""
    for line in file.readlines():
        if line[0] == '>':
            sequence_labels.append(line[1:].replace('\n', ''))
            if(seq != ""):
                sequences.append([int(s) for s in seq])
                seq = ""
        else:
            seq += line.replace('\n', '')
    sequences.append([int(s) for s in seq])
    file.close()

    num_sequences = len(sequences)
    num_cols = len(sequences[0])

    print("seqs: ", num_sequences, " cols: ", num_cols)

    genes = np.zeros((num_sequences, num_cols), dtype=int)
    for i in range(num_sequences):
        for j in range(num_cols):
            if(sequences[i][j] == '-'):
                genes[i][j] = -1
            else:
                genes[i][j] = int(sequences[i][j])

def get_subsets(M):
    M_copy = M[:]

def powerset(M, w, s):
    power_set = []

    for top_index in range(len(M)):
        low_bound = max(0, top_index - w + 1)
        c = [list(subset) + [M[top_index]] for r in range(1, min(s, top_index-low_bound + 1)) for subset in list(combinations(M[low_bound:top_index], r))]
        power_set.extend(c)
    
    return power_set

def haplotype_bound(M):
    haplotypes = np.unique(M, axis = 0)

    return len(haplotypes) - M.shape[1] - 1 # Assume that input cleaned so columns are all nonzero

def RecMin(M, w, s):
    power_set = powerset(range(M.shape[1]), w, s)
    recombs = []

    for subset in power_set:
        # recall that it is sorted by largest element in list
        b = haplotype_bound(M[:, subset])

        # should have b recombinations within interval spanning the subset.
        lowest = np.min(subset)
        highest = np.max(subset)
        current_in_range = len([recomb for recomb in recombs if recomb > lowest])

        if current_in_range < b:
            # Need to add more recombs. A recombination at x means just before site x
            recombs.extend([highest for i in range(b - current_in_range)])

    return len(recombs)



def main(argv):
    inputfile = ''
    w = 0
    s = 0

    try:
        opts, args = getopt.getopt(argv, "hi:w:s:", ["ifile="])
    except getopt.GetoptError:
        print('haplotype_bound.py -i <input file> -w <width> -s <size>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('haplotype_bound.py -i <input file> -w <width> -s <size>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-w"):
            w = int(arg)
        elif opt in ("-s"):
            s = int(arg)

    print('Input file is: ', inputfile)
    print('w: ', w, ' s:, ', s)

    read_sequences(inputfile)

    print("H: ", haplotype_bound(genes))
    print("H+: ", RecMin(genes, w, s))


if __name__ == "__main__":
    main(sys.argv[1:])
