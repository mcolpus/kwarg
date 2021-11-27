import numpy as np



do_clean = True


sequences = []

filename = "mafft alligned - hcov_global_about_500.txt"
file = open(filename, "r")
seq = ""
for line in file.readlines():
    if line[0] == '>':
        if(seq != ""):
            sequences.append(seq)
            seq = ""
    else:
        seq += line.replace('\n', '')
file.close()

num_sequences = len(sequences)
num_cols = len(sequences[0])

print("seqs: ", num_sequences)
print("cols: ", num_cols)

# Convert to binary sequences. 0 will be most common, 1 for everything else, -1 if - or unknown
# Will have each site as a row then transpose
genes = []

bases = ['a', 'g', 'c', 't']

col_out = 0
for col_in in range(num_cols):
    counts = {}
    counts['a'] = 0
    counts['t'] = 0
    counts['c'] = 0
    counts['g'] = 0
    for j in range(num_sequences):
        c = sequences[j][col_in]
        if(c in bases):
            counts[c] += 1

    max_key = max(counts, key=counts.get)
    
    if do_clean:
        counts[max_key] = 0
        if(max(counts.values()) <= 1):
            # all mutations in this column are unique so can remove
            continue

    site = [0] * num_sequences
    for j in range(num_sequences):
        c = sequences[j][col_in]
        if(c == max_key):
            site[j] = 0
        elif(c in bases):
            site[j] = 1
        else:
            site[j] = -1
    genes.append(site)

genes = np.array(genes, dtype=int)
genes = genes.transpose()
print(genes)

print(genes.shape)

# clean: 




# output to txt

output = "binary_sequences"

file = open(output, 'w')
for j in range(genes.shape[0]):
    line = ""
    for i in range(genes.shape[1]):
        if(genes[j][i] == -1):
            line += '-'
        else:
            line += str(genes[j][i])
    file.write(line + "\n")

file.close()