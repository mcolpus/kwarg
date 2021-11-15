import numpy as np

# Read in input file line at a time to give each of the gene sequences

sequences = [] # ["00", "10", "01", "11"]

# sequences = ["000000000", "000111000", "111000111", "111111111"]

filename = "kreitman_snp.txt"
file1 = open(filename, "r") 
for line in file1.readlines():
    sequences.append(line.replace('\n', ''))
file1.close()

num_sequences = len(sequences)
num_cols = len(sequences[0])


for s in sequences:
    print(s)

genes = np.zeros((num_sequences, num_cols), dtype=int)
for i in range(num_sequences):
    for j in range(num_cols):
        genes[i][j] = int(sequences[i][j])




def incomp(i, j):
    # want to check is col i and col j are incompatible
    b00 = False
    b01 = False
    b10 = False
    b11 = False

    for k in range(num_sequences):
        a = genes[k][i]
        b = genes[k][j]
        if(a == 0 and b == 0):
            b00 = True
        if(a == 1 and b == 0):
            b10 = True
        if(a == 0 and b == 1):
            b01 = True
        if(a == 1 and b == 1):
            b11 = True
    
    return (b00 and b01 and b10 and b11)

def hudson():
    regions = []
    i = 0
    j = 1
    while j < num_cols:
        for a in range(i, j):
            if incomp(a, j):
                regions.append((a,j))
                i = j
                break
        j += 1
    return (len(regions), regions)


(recombs, regions) = hudson()

print("hudson recombs: ", recombs)
for (a, b) in regions:
    print("({}, {})".format(a, b))


def two_hudson():
    # will have a and b incomp with i and j
    regions = []
    left = 0
    i = 2
    while i < num_cols - 1:
        ab_incomps = []
        for a in range(left, i):
            if(incomp(a, i)):
                ab_incomps.append(a)
        if(len(ab_incomps) >= 2):
            for j in range(i+1, num_cols):
                x = 0
                a = -1
                b = -1
                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], j)):
                        if(a == -1):
                            a = ab_incomps[x]
                        elif(b == -1):
                            b = ab_incomps[x]
                            break
                        else:
                            print("Error, shouldn't get here")
                    x += 1
                
                if(a != -1 and b != -1):
                    regions.append((a,b,i,j))
                    left = i
                    break
        
        i += 1

    return (len(regions), regions)

def find_incomp_subset(set, target):
    subset = []
    for i in set:
        if incomp(i, target):
            subset.append(i)
    return subset

def find_b_r(k, r, b_index, current_incomps):
    # have found b1..b_(r-1) (b_r = b_index) and know current_incomps are incompatible with all of them
    for br in range(b_index, num_cols):
        subset = find_incomp_subset(current_incomps, br)
        if(len(subset) >= k):
            if(r == k):
                # done
                return ([b_index], subset[0:k])
            else:
                (bs, ays) = find_b_r(k, r+1, br + 1, subset)
                if(len(bs) > 0):
                    bs.append(br)
                    return(bs, ays)

    return ([], [])

print(incomp(3, 6))       

def k_hudson(k):
    # will have a1 to ak incompatible with b1 to bk
    regions = []
    left = 0
    b1 = k
    while b1 <= num_cols - k:
        print("b1: ", b1, " left: ", left)
        a_b1_incomps = []
        for a in range(left, b1):
            if(incomp(a, b1)):
                a_b1_incomps.append(a)
        if(len(a_b1_incomps) >= k):
            # Now need to find rest of b2...bk which will work
            # Use a recursive function to find next b (and all the rest)
            (b_list, a_list) = find_b_r(k, 2, b1+1, a_b1_incomps)
            if(len(b_list) > 0):
                b_list.append(b1)
                b_list.reverse()
                regions.append((a_list, b_list))
                left = b1
                b1 += k-1
        
        b1 += 1

    return (len(regions), regions)

(recombs, regions) = two_hudson()

print("2-hudson recombs: ", recombs)
for (a, b, i , j) in regions:
    print("({}, {}, {}, {})".format(a, b, i, j))

(recombs, regions) = k_hudson(3)

print("3-hudson recombs: ", recombs)
for (ays, bs) in regions:
    print("ays: ", ays)
    print("bs: ", bs)




