import numpy as np
from numpy.lib.function_base import diff

sequences = []

filename = "binary_sequences.txt"
file = open(filename, "r")
seq = ""
for line in file.readlines():
    sequences.append(line.replace('\n', '').replace('-', '9'))
file.close()

sequences = sequences[:300]

num_seqs = len(sequences)
num_cols = len(sequences[0])
print("seqs: ", num_seqs)
print("cols: ", num_cols)

seq_array = np.zeros((num_seqs, num_cols), dtype=int)
for i in range(num_seqs):
    for j in range(num_cols):
        seq_array[i][j] = int(sequences[i][j])

print("have converted to np array")
distances = np.zeros((num_seqs, num_seqs), dtype=int)
diffs_data = []

for i in range(num_seqs):
    if i % 10 == 0:
        print(i)
    for j in range(i+1, num_seqs):
        diffs = 0
        for c in range(num_cols):
            if sequences[i][c] != sequences[j][c]:
                diffs += 1
        distances[i][j] = diffs
        distances[j][i] = diffs
        diffs_data.append(diffs)


# print(diffs_data)
print("avg", sum(diffs_data) / len(diffs_data))
print("max", max(diffs_data))

ones = [0] * num_seqs
for i in range(num_seqs):
    for c in range(num_cols):
        if sequences[i][c] == '1':
            ones[i] += 1
print(ones)
