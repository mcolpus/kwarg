from typing import Sequence
import numpy as np


do_clean = True
num_states = 2
keep_seq_labels = True
use_ref = True  # Reference should be first sequence


sequences = []
sequence_labels = []

filename = "eng_alignment_cleaned.fasta"
file = open(filename, "r")
seq = ""
for line in file.readlines():
    if line[0] == '>':
        sequence_labels.append(line[1:].replace('\n', ''))
        if(seq != ""):
            sequences.append(seq.upper())
            seq = ""
    else:
        seq += line.replace('\n', '')
sequences.append(seq)
file.close()

num_sequences = len(sequences)
num_cols = len(sequences[0])

print("seqs: ", num_sequences)
print("cols: ", num_cols)

bases = ['A', 'G', 'C', 'T']
opposite = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}
not_opposite = {  # This is just used when turning into 4 states with a reference
    'A': 'G',
    'G': 'A',
    'T': 'C',
    'C': 'T'
}
bin_sequences = []

if use_ref:
    ref_seq = sequences[0]
    ref_label = sequence_labels[0]
    # sequences.remove(ref_seq)
    # sequence_labels.remove(ref_label)

    if num_states == 2:
        for seq in sequences:
            bin_seq = [0] * num_cols
            for i in range(num_cols):
                if(seq[i] == ref_seq[i]):
                    bin_seq[i] = 0
                elif(seq[i] in bases):
                    bin_seq[i] = 1
                else:
                    bin_seq[i] = 9  # Means error
            bin_sequences.append(bin_seq)
    elif num_states == 3:
        for seq in sequences:
            bin_seq = [0] * num_cols
            for i in range(num_cols):
                if seq[i] == ref_seq[i]:
                    bin_seq[i] = 0
                elif seq[i] == opposite[ref_seq[i]]:
                    bin_seq[i] = 1
                elif seq[i] in bases:
                    bin_seq[i] = 2
                else:
                    bin_seq[i] = 9  # Means error
            bin_sequences.append(bin_seq)
    elif num_states == 4:
        for seq in sequences:
            bin_seq = [0] * num_cols
            for i in range(num_cols):
                if seq[i] == ref_seq[i]:
                    bin_seq[i] = 0
                elif seq[i] == opposite[ref_seq[i]]:
                    bin_seq[i] = 1
                elif seq[i] == not_opposite[ref_seq[i]]:
                    bin_seq[i] = 2
                elif seq[i] in bases:
                    bin_seq[i] = 3
                else:
                    bin_seq[i] = 9  # Means error
            bin_sequences.append(bin_seq)

    bin_sequences = np.array(bin_sequences, dtype=int)

else:
    # Convert to binary sequences. 0 will be most common, 1 for everything else, 9 if - or unknown
    # Will have each site as a row then transpose
    bin_cols = []

    for col in range(num_cols):
        counts = {}
        counts['A'] = 0
        counts['T'] = 0
        counts['C'] = 0
        counts['G'] = 0
        for j in range(num_sequences):
            c = sequences[j][col]
            if(c in bases):
                counts[c] += 1

        max_key = max(counts, key=counts.get)

        site = [0] * num_sequences
        for j in range(num_sequences):
            c = sequences[j][col]

            if num_states == 2:
                if(c == max_key):
                    site[j] = 0
                elif(c in bases):
                    site[j] = 1
                else:
                    site[j] = 9
            elif num_states == 3:
                if(c == max_key):
                    site[j] = 0
                elif(c == opposite[max_key]):
                    site[j] = 1
                elif(c in bases):
                    site[j] = 2
                else:
                    site[j] = 9
            elif num_states == 4:
                if(c == max_key):
                    site[j] = 0
                elif(c == opposite[max_key]):
                    site[j] = 1
                elif(c == not_opposite[max_key]):
                    site[j] = 2
                elif(c in bases):
                    site[j] = 3
                else:
                    site[j] = 9
        bin_cols.append(site)

    bin_cols = np.array(bin_cols, dtype=int)
    bin_sequences = bin_cols.transpose()

print(bin_sequences)

print(bin_sequences.shape)

# clean:
print("Finished reading in data, now cleaning")
if do_clean:
    bin_cols = bin_sequences.transpose()
    empty_cols = []
    for i in range(len(bin_cols)):
        others_seen = 0
        for c in bin_cols[i]:
            if c != 0 and c != 9:
                others_seen += 1
                if others_seen >= 2:
                    break
        if others_seen < 2:
            empty_cols.append(i)
    
    # Now delete all these columns
    bin_cols = np.delete(bin_cols, empty_cols, axis=0)
    bin_sequences = bin_cols.transpose()

# output to txt
print("Finished cleaning, now outputting")
output = "binary_sequences"

file = open(output, 'w')
for j in range(bin_sequences.shape[0]):
    if keep_seq_labels:
        file.write(">" + sequence_labels[j] + "\n")
    line = ""
    for i in range(bin_sequences.shape[1]):
        if(bin_sequences[j][i] == 9):
            line += '-'
        else:
            line += str(bin_sequences[j][i])
    file.write(line + "\n")

file.close()
