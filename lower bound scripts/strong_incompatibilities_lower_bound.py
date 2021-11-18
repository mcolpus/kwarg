import numpy as np

# Read in input file line at a time to give each of the gene sequences

print_regions = True

run_hudson = False
run_2_hudson = False
run_3_hudson = False

sequences = [] # ["00", "10", "01", "11"]

# sequences = ["000000000", "000111000", "111000111", "111111111"]

filename = "binary_sequences"
file1 = open(filename, "r") 
for line in file1.readlines():
    sequences.append(line.replace('\n', ''))
file1.close()

num_sequences = len(sequences)
num_cols = len(sequences[0])


# for s in sequences[0:2]:
#     print(s)

genes = np.zeros((num_sequences, num_cols), dtype=int)
for i in range(num_sequences):
    for j in range(num_cols):
        if(sequences[i][j] == '-'):
            genes[i][j] = -1
        else:
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
        if(not b00 and a == 0 and b == 0):
            b00 = True
            if b00 and b01 and b10 and b11:
                break
        if(not b10 and a == 1 and b == 0):
            b10 = True
            if b00 and b01 and b10 and b11:
                break
        if(not b01 and a == 0 and b == 1):
            b01 = True
            if b00 and b01 and b10 and b11:
                break
        if(not b11 and a == 1 and b == 1):
            b11 = True
            if b00 and b01 and b10 and b11:
                break
    
    return b00 and b01 and b10 and b11

def incomp_region(group_1, group_2):
    for i in group_1:
        for j in group_2:
            if(not incomp(i,j)):
                print(i, "and ", j, " are not incompatible")

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
                return ([br], subset[0:k])
            else:
                (bs, ays) = find_b_r(k, r+1, br + 1, subset)
                if(len(bs) > 0):
                    bs.append(br)
                    return(bs, ays)

    return ([], [])

def k_hudson(k):
    # will have a1 to ak incompatible with b1 to bk
    regions = []
    left = 0
    b1 = k
    while b1 <= num_cols - k:
        a_b1_incomps = []
        for a in range(left, b1):
            if(incomp(a, b1)):
                a_b1_incomps.append(a)
        if(len(a_b1_incomps) >= k):
            print("b1: ", b1, " incomps: ", a_b1_incomps)
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

def intersect(l1, l2):
    return [i for i in l1 if i in l2]

_far_left_boundary = 0
_close_left_boundary = 0
_incomps_with_earlier = []

def k_hudson_non_overlap_iterate(k, r, b_prev, b_k, shared_incomps):
    # print(" iter r {}  b_prev {}  b_k {}  left {}  shared {}".format(r, b_prev, b_k, _close_left_boundary, shared_incomps))
    global _close_left_boundary
    global _incomps_with_earlier

    if(r==k):
        # print("r==k shared: ", shared_incomps[0:k-1] + shared_incomps[-1:])
        return (shared_incomps[0:k-1] + shared_incomps[-1:], [b_k])
    else:
        for b_r in range(b_prev + 1, b_k):
            intersection = intersect(shared_incomps, _incomps_with_earlier[b_r])
            if len(intersection) >= k and max(intersection) >= _close_left_boundary:
                (a_list, b_list) = k_hudson_non_overlap_iterate(k, r+1, b_r, b_k, intersection)
                if len(b_list) > 0:
                    b_list.append(b_r)
                    return (a_list, b_list)

    return ([], [])



def k_hudson_non_overlap(k):
    # will have a1 to ak incompatible with b1 to bk
    # Given a's incomp with b's, and x's with y's
    # Then must have a's...b1 - b&x's - xk .. y's
    # This is much harder to be optimal! So will try to minimize b_k, and then b1
    regions = []
    global _far_left_boundary
    _far_left_boundary = 0 # a1...a(k-1) need to be equal or right of this
    global _close_left_boundary
    _close_left_boundary = k-1 # ak needs to be equal or right of this!
    global _incomps_with_earlier
    _incomps_with_earlier = [[] for b in range(num_cols)] 

    b = k
    while b < num_cols:
        
        ealier_incomps = []
        for a in range(_far_left_boundary, b):
            if incomp(a, b):
                ealier_incomps.append(a)
        _incomps_with_earlier[b] = ealier_incomps
        # print("b: ", b, " far-left: ", far_left_boundary, " incomp with ", ealier_incomps)

        if len(ealier_incomps) >= k and max(ealier_incomps) >= _close_left_boundary:
            (a_list, b_list) = k_hudson_non_overlap_iterate(k, 1, _close_left_boundary, b, ealier_incomps)
            if len(b_list) > 0:
                b_list.reverse()
                regions.append((a_list, b_list))
                print("found: ", a_list, "  ", b_list)

                _far_left_boundary = b_list[0]
                _close_left_boundary = b_list[-1]
                b = _close_left_boundary
                continue

        b += 1

    return (len(regions), regions)


if run_hudson:
    (recombs, regions) = hudson()

    print("hudson recombs: ", recombs)

    if print_regions:
        for (a, b) in regions:
            print("({}, {})".format(a, b))

if run_2_hudson:
    (recombs, regions) = two_hudson()

    print("2-hudson recombs: ", recombs)

    if print_regions:
        for (a, b, i , j) in regions:
            print("({}, {}, {}, {})".format(a, b, i, j))

if run_3_hudson:
    (recombs, regions) = k_hudson(3)

    print("3-hudson recombs: ", recombs)

    if print_regions:
        for (ays, bs) in regions:
            print("ays: ", ays, "   bs: ", bs)


(recombs, regions) = k_hudson_non_overlap(2)

print("2-hudson improved recombs: ", recombs)

if print_regions:
    for (ays, bs) in regions:
        print("ays: ", ays, "   bs: ", bs)


(recombs, regions) = k_hudson_non_overlap(3)

print("3-hudson improved recombs: ", recombs)

if print_regions:
    for (ays, bs) in regions:
        print("ays: ", ays, "   bs: ", bs)




