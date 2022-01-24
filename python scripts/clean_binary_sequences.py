
import enum
import sys
import getopt

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

    num_sequences = len(sequences)
    num_cols = len(sequences[0])

    print("seqs: ", num_sequences)
    print("cols: ", num_cols)

    return (sequence_labels, sequences)

def remove_uninformative_cols(sequences):
    seq_length = len(sequences[0])
    keep_columns = []

    for site in range(seq_length):
        number_of_ones = 0
        for seq in sequences:
            if seq[site] == 1:
                number_of_ones += 1
                if number_of_ones >= 2:
                    keep_columns.append(site)
                    break
    
    print("remove uniformative keeps", keep_columns)
    
    if len(keep_columns) == seq_length:
        return (False, sequences)
    else:
        return (True, [[seq[site] for site in keep_columns] for seq in sequences])

def merge_adj_identical_cols(sequences):
    seq_length = len(sequences[0])
    if seq_length == 0:
        return (False, sequences)

    keep_columns = [0] # Keep columns which are not identical to the previous

    for site in range(1, seq_length):
        for seq in sequences:
            if seq[site] != seq[site-1]:
                keep_columns.append(site)
                break
    

        
    print("merge_adj_identical_cols keeps", keep_columns)
    
    if len(keep_columns) == seq_length:
        return (False, sequences)
    else:
        return (True, [[seq[site] for site in keep_columns] for seq in sequences])

def remove_identical_rows(sequences, sequence_labels):
    rows_to_be_kept = []

    row_descriptions = {}

    for i, seq in enumerate(sequences):
        description = ""
        for site, v in enumerate(seq):
            if v == 1:
                description += str(site) + ","
        
        if description not in row_descriptions:
            row_descriptions[description] = 1
            rows_to_be_kept.append(i)
    
    print("remove_identical_rows keeping", rows_to_be_kept)

    if len(rows_to_be_kept) == len(sequences):
        return (False, sequences, sequence_labels)
    else:
        new_seqs = [sequences[i] for i in rows_to_be_kept]
        new_labels = [sequence_labels[i] for i in rows_to_be_kept]

        return (True, new_seqs, new_labels)

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

    (sequence_labels, sequences) = read_sequences(inputfile)

    somethings_changed = True
    while(somethings_changed):
        print("run iteration")
        (change1, sequences) = remove_uninformative_cols(sequences)
        (change2, sequences) = merge_adj_identical_cols(sequences)
        (change3, sequences, sequence_labels) = remove_identical_rows(sequences, sequence_labels)
        somethings_changed = change1 or change2 or change3
    
    print("done")
    output(sequences, sequence_labels, outputfile)



if __name__ == "__main__":
    main(sys.argv[1:])