
import sys
import getopt
import random

verbose = False

def read_sequences(inputfile):
    sequences = []
    sequence_labels = []

    file = open(inputfile, "r")
    seq = ""
    first_loop = True
    has_labels = True

    for line in file.readlines():
        if first_loop:
            first_loop = False
            has_labels = line[0] == '>'
        
        if has_labels:
            if line[0] == '>':
                sequence_labels.append(line[1:].replace('\n', ''))
                if(seq != ""):
                    sequences.append([int(s) for s in seq])
                    seq = ""
            else:
                seq += line.replace('\n', '')
        else:
            sequences.append([int(s) for s in line.replace('\n', '')])
    
    if has_labels:
        sequences.append([int(s) for s in seq])
    file.close()

    if len(sequences) == 0:
        print("No sequences in input!")
        return ([], [], 0)
    
    if not has_labels:
        sequence_labels = ["seq {}".format(i) for i in range(len(sequences))]
    num_cols = len(sequences[0])

    return (sequence_labels, sequences, num_cols)

def split_random(sequence_labels, sequences, sample_size):
    num_sequences = len(sequences)

    indexes = [i for i in range(num_sequences)]
    random.shuffle(indexes)

    samples_indexes = list()


    for i in range(0, num_sequences//sample_size):
        samples_indexes.append(indexes[i*sample_size:i*sample_size + sample_size])
    
    print(samples_indexes)

    samples = list()
    for sample_indexes in samples_indexes:
        sample_seq = list()
        sample_labels = list()
        for i in sample_indexes:
            sample_seq.append(sequences[i])
            sample_labels.append(sequence_labels[i])
        samples.append((sample_labels, sample_seq))

        print(sample_labels)

    return samples


def output(binary_sequences, sequence_labels, outfile):
    file = open(outfile, 'w')
    for j in range(len(binary_sequences)):
        file.write(">" + sequence_labels[j] + "\n")
        line = ""
        for v in binary_sequences[j]:
            line += str(v)
        file.write(line + "\n")

    file.close()

def main(argv):
    inputfile = ''
    outputfile = ''
    sample_size = 10

    try:
        opts, args = getopt.getopt(argv, "hi:o:s:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('clean_binary_sequences.py -i <input file> -o <output file> -s <size to split into>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('clean_binary_sequences.py -i <input file> -o <output file> -s <size to split into>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-s"):
            sample_size = int(arg)

    print('Input file is: ', inputfile)
    print('Output file is: ', outputfile)
    print('Sample size: ', sample_size)

    (sequence_labels, sequences, num_cols) = read_sequences(inputfile)

    num_sequences = len(sequences)
    print("has ", num_sequences, " sequences")
    
    samples = split_random(sequence_labels, sequences, sample_size)

    for i, (sample_labels, sample_seq) in enumerate(samples):
        output(sample_seq, sample_labels, outputfile + str(i) + ".txt")



if __name__ == "__main__":
    main(sys.argv[1:])