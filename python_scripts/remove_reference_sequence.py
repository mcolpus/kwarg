import sys
import getopt

verbose = False

def read_sequences(inputfile):
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

    return (sequence_labels, sequences)


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

    try:
        opts, args = getopt.getopt(argv, "hm:i:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('remove_reference_sequence.py -i <input file> -o <output file>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('remove_reference_sequence.py -i <input file> -o <output file>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print('Input file is: ', inputfile)
    print('Output file is: ', outputfile)

    (sequence_labels, sequences) = read_sequences(inputfile)
    output(sequences[1:], sequence_labels[1:], outputfile)


if __name__ == "__main__":
    main(sys.argv[1:])