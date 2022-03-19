from re import sub
from xmlrpc.client import Unmarshaller
import numpy as np
import sys
import getopt
import random

from clean_binary_sequences import read_sequences

genes = np.zeros(0)
incomp_grid = np.zeros(0)
num_sequences = 0
num_cols = 0


def read_genes(inputfile):
    global genes
    global num_sequences
    global num_cols

    (_, sequences, num_cols) = read_sequences(inputfile)
    num_sequences = len(sequences)

    print("seqs: ", num_sequences, " cols: ", num_cols)

    genes = np.zeros((num_sequences, num_cols), dtype=int)
    for i in range(num_sequences):
        for j in range(num_cols):
            if(sequences[i][j] == '-'):
                genes[i][j] = -1
            else:
                genes[i][j] = int(sequences[i][j])


def incomp(i, j):
    global genes
    global num_sequences
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


def make_incomp_grid():
    global num_cols
    global incomp_grid

    incomp_grid = np.zeros((num_cols, num_cols))
    for i in range(num_cols):
        for j in range(num_cols):
            incomp_grid[i][j] = incomp(i, j)


def hudson(overlap=True):
    regions = []
    i = 0
    j = 1
    while j < num_cols:
        for a in range(i, j):
            if incomp(a, j):
                regions.append((a, j))
                i = j
                if not overlap:
                    i += 1
                    j += 1
                break
        j += 1
    return (len(regions), regions)


def two_hudson():
    # will have a and b incomp with i and j.
    # Can have some overlap: a1-b1=i1 a2 j1 b2=i2-j2, but can't allow overlap of numbers
    # Idea is that to remove any of these recombinations needs two recurrent mutations
    regions = []
    left = 0
    i = 2
    prev_j = -1
    while i < num_cols - 1:
        ab_incomps = []
        for a in range(left, i):
            if a == prev_j:
                continue  # Can't have a or b same as previous j
            if(incomp(a, i)):
                ab_incomps.append(a)

        if(len(ab_incomps) >= 2):
            for j in range(i+1, num_cols):
                x = 0
                a = -1
                b = -1

                # First try an find a, this is allowed to be before previous j
                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], j)):
                        a = ab_incomps[x]
                        x += 1
                        break
                    x += 1

                # Now look for b which must be after previous j
                while (x < len(ab_incomps) and ab_incomps[x] < prev_j):
                    x += 1

                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], j)):
                        b = ab_incomps[x]
                        break
                    x += 1

                if(a != -1 and b != -1):
                    regions.append((a, b, i, j))
                    left = i + 1
                    i = j + 1  # i will increase once more at end of iteration
                    prev_j = j
                    break

        i += 1

    return (len(regions), regions)


def two_hudson2():
    # This version aims to make j minimal rather than i
    regions = []
    left = 0
    j = 3
    prev_j = -1
    while j < num_cols:
        ab_incomps = [a for a in range(
            left, j-1) if a != prev_j and incomp(a, j)]

        if(len(ab_incomps) >= 2):
            for i in range(j-1, max(prev_j + 2, ab_incomps[1]+1), -1):
                x = 0
                a = -1
                b = -1

                # First try an find a, this is allowed to be before previous j
                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], i)):
                        a = ab_incomps[x]
                        x += 1
                        break
                    x += 1

                # Now look for b which must be after previous j
                while (x < len(ab_incomps) and ab_incomps[x] < prev_j):
                    x += 1

                while x < len(ab_incomps):
                    if(incomp(ab_incomps[x], i)):
                        b = ab_incomps[x]
                        break
                    x += 1

                if(a != -1 and b != -1):
                    regions.append((a, b, i, j))
                    left = i + 1
                    prev_j = j
                    j += 2  # j will increase once more at end of iteration
                    break

        j += 1

    return (len(regions), regions)


def find_br_down_no_overlap(k, r, start_search, current_incomps):
    # have found bk..b_(r+1) (b_r = b_index) and know current_incomps are incompatible with all of them

    for br in range(start_search, 0, -1):
        subset = [i for i in current_incomps if incomp(i, br)]
        if(len(subset) >= k and subset[k-1] <= br - r):
            if(r == 1):
                # done
                return ([br], subset[0:k])
            else:
                (bs, ays) = find_br_down_no_overlap(k, r-1, br - 1, subset)
                if(len(bs) > 0):
                    bs.append(br)
                    return(bs, ays)

    return ([], [])


def k_hudson_no_overlap(k):
    # will have a1 to ak incompatible with b1 to bk.
    # Each region must be completely separate
    # To make them smallest we actually base it on bk
    regions = []
    left = 0
    bk = 2*k - 1
    while bk < num_cols:
        a_incomps = [a for a in range(left, bk) if incomp(a, bk)]
        if(len(a_incomps) >= k and a_incomps[k-1] <= bk - k):
            # Checks if there is enough of a gap to fit b1,...,b(k-1)

            (b_list, a_list) = find_br_down_no_overlap(k, k-1, bk-1, a_incomps)
            if(len(b_list) > 0):
                b_list.append(bk)
                regions.append((a_list, b_list))
                left = bk+1
                bk = left + 2*k - 2

        bk += 1

    return (len(regions), regions)


# Tail overlap section is problematic as it allows too much

def intersect(l1, l2):
    return [i for i in l1 if i in l2]


_far_left_boundary = 0
_close_left_boundary = 0
_incomps_with_earlier = []


def k_hudson_tail_overlap_iterate(k, r, b_prev, b_k, shared_incomps):
    # print(" iter r {}  b_prev {}  b_k {}  left {}  shared {}".format(r, b_prev, b_k, _close_left_boundary, shared_incomps))
    global _close_left_boundary
    global _incomps_with_earlier

    if(r == k):
        # print("r==k shared: ", shared_incomps[0:k-1] + shared_incomps[-1:])
        return (shared_incomps[0:k-1] + shared_incomps[-1:], [b_k])
    else:
        for b_r in range(b_prev + 1, b_k):
            intersection = intersect(
                shared_incomps, _incomps_with_earlier[b_r])
            if len(intersection) >= k and max(intersection) >= _close_left_boundary:
                (a_list, b_list) = k_hudson_tail_overlap_iterate(
                    k, r+1, b_r, b_k, intersection)
                if len(b_list) > 0:
                    b_list.append(b_r)
                    return (a_list, b_list)

    return ([], [])


def k_hudson_tail_overlap(k):
    # will have a1 to ak incompatible with b1 to bk
    # Given a's incomp with b's, and x's with y's
    # Then must have a's...b1 - b&x's - xk .. y's
    # This is much harder to be optimal! So will try to minimize b_k, and then b1
    regions = []
    global _far_left_boundary
    _far_left_boundary = 0  # a1...a(k-1) need to be equal or right of this
    global _close_left_boundary
    _close_left_boundary = k-1  # ak needs to be equal or right of this!
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
            (a_list, b_list) = k_hudson_tail_overlap_iterate(
                k, 1, _close_left_boundary, b, ealier_incomps)
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


def masked_hudson(mask):
    global num_cols
    global incomp_grid

    unmasked = [a for a in range(num_cols) if not (a in mask)]

    regions = []
    left = 0

    for j in unmasked:
        for i in range(left, j):
            if i in mask:
                continue
            if incomp_grid[i][j]:
                regions.append((i, j))
                left = j
                break

    return (len(regions), regions)


def split_out_dups(sites_used):
    sites_used_no_dups = []
    sites_used_dups = []
    # how many site masking events could remove two recombinations
    potential_dups_removable = 0

    prev_was_removable_dup = False
    i = 1 
    while i < len(sites_used):
        if sites_used[i-1] == sites_used[i]:
            sites_used_dups.append(sites_used[i])
            i += 2
            if not prev_was_removable_dup:
                potential_dups_removable += 1
                prev_was_removable_dup = True
            else:
                prev_was_removable_dup = False
        else:
            sites_used_no_dups.append(sites_used[i-1])
            prev_was_removable_dup = False
            i += 1
    sites_used_no_dups.append(sites_used[-1])

    print("no_dups: ", sites_used_no_dups)
    print("dups: ", sites_used_dups)
    print("removable: ", potential_dups_removable)


_target_ratio = 2.0

def crude_adverserial_hudson_iter(current_mask, m, ignored_sites):
    global _target_ratio
    # can place m more recurrent mutations
    (num_regions, regions) = masked_hudson(current_mask)
    if (m == 0 or num_regions == 0):
        return (num_regions, regions, current_mask)

    sites_used = [i for ij in regions for i in ij]
    # may have duplicates but not a problem

    best_num_regions = sys.maxsize
    best_regions = []
    best_mask = []

    target = max(0, float(num_regions) - (2 * m)) # Best that can be done is to 
    sites_considered = 0
    for site in sites_used:
        if site in ignored_sites:
            continue

        (new_num_regions, new_regions, new_mask) = crude_adverserial_hudson_iter(current_mask + [site], m-1, ignored_sites[:])
        if (new_num_regions < best_num_regions):
            # Adverserial so want to minimize
            best_num_regions = new_num_regions
            best_regions = new_regions
            best_mask = new_mask

            if best_num_regions <= target:
                break

        ignored_sites.append(site)

    return [best_num_regions, best_regions, best_mask]


def crude_adverserial_hudson(m):
    make_incomp_grid()
    # adverserial is allowed to select m sites to place recurrent mutations on
    return crude_adverserial_hudson_iter([], m, [])

def target_adverserial_hudson_iter(current_mask, m, ignored_sites):
    global _target_ratio
    # can place m more recurrent mutations
    (num_regions, regions) = masked_hudson(current_mask)
    if (m == 0 or num_regions == 0):
        return (num_regions, regions, current_mask)

    sites_used = [i for ij in regions for i in ij]
    # may have duplicates but not a problem

    sites_used_no_dups = []
    sites_used_dups = []
    # how many site masking events could remove two recombinations
    potential_dups_removable = 0

    prev_was_removable_dup = False
    i = 1
    while i < len(sites_used):
        if sites_used[i-1] == sites_used[i]:
            sites_used_dups.append(sites_used[i])
            i += 2
            if not prev_was_removable_dup:
                potential_dups_removable += 1
                prev_was_removable_dup = True
            else:
                prev_was_removable_dup = False
        else:
            sites_used_no_dups.append(sites_used[i-1])
            prev_was_removable_dup = False
            i += 1
    sites_used_no_dups.append(sites_used[-1])

    best_num_regions = sys.maxsize
    best_regions = []
    best_mask = []

    target = num_regions - (2 * min(m, potential_dups_removable)) - max(0, m - potential_dups_removable)
    sites_considered = 0
    for site in sites_used_dups + sites_used_no_dups:
        if site in ignored_sites:
            continue

        (new_num_regions, new_regions, new_mask) = target_adverserial_hudson_iter(current_mask + [site], m-1, ignored_sites[:])
        if (new_num_regions < best_num_regions):
            # Adverserial so want to minimize
            best_num_regions = new_num_regions
            best_regions = new_regions
            best_mask = new_mask

            if best_num_regions <= target:
                break

        ignored_sites.append(site)

    return [best_num_regions, best_regions, best_mask]


def target_adverserial_hudson(m):
    make_incomp_grid()
    # adverserial is allowed to select m sites to place recurrent mutations on
    return target_adverserial_hudson_iter([], m, [])


def random_adverserial_hudson_iter(current_mask, m):
    # can place m more recurrent mutations
    (num_regions, regions) = masked_hudson(current_mask)
    if (m == 0 or num_regions == 0):
        return (num_regions, regions, current_mask)

    sites_used = [i for ij in regions for i in ij]

    site = random.choice(sites_used)

    return random_adverserial_hudson_iter(current_mask + [site], m-1)


def random_adverserial_hudson(m):
    # adverserial is allowed to select m sites to place recurrent mutations on
    make_incomp_grid()

    best_num_regions = sys.maxsize
    best_regions = []
    bset_mask = []

    for i in range(200):
        (num_regions, regions, mask) = random_adverserial_hudson_iter([], m)
        if (num_regions < best_num_regions):
            best_num_regions = num_regions
            best_regions = regions
            best_mask = mask

    return [best_num_regions, best_regions, best_mask]


def main(argv):
    inputfile = ''
    ks = [1]
    print_regions = False
    overlap = True
    adverserial = False

    try:
        opts, args = getopt.getopt(argv, "hi:k:pna", ["ifile="])
    except getopt.GetoptError:
        print('hudson_bound.py -i <input file> -k <k-incompatibilities or m rec-muts> -p (print regions) -n (non-overlap) -a (adversial approach)')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(
                'hudson_bound.py -i <input file> -k <find k-incompatibilities> -p (print regions) -n (non-overlap)')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-k"):
            ks = [int(a) for a in arg.split(",")]
        elif opt in ("-p"):
            print_regions = True
        elif opt in ("-n"):
            overlap = False
        elif opt in ("-a"):
            adverserial = True

    print('Input file is: ', inputfile)
    print('k2 is: ', ks)

    read_genes(inputfile)
    # Input gets read to global variables

    for k in ks:
        if adverserial:
            (recombs, regions, mask) = crude_adverserial_hudson(k)
            print("crude adverserial recombs: ", recombs,
                  " with m: ", k, " and mask: ", sorted(mask))
            print(regions)

            (recombs, regions, mask) = target_adverserial_hudson(k)
            print("target adverserial recombs: ", recombs,
                  " with m: ", k, " and mask: ", sorted(mask))
            print(regions)

            if print_regions:
                for (a, b) in regions:
                    print("({}, {})".format(a, b))
            continue

        if (k == 1):
            (recombs, regions) = hudson(overlap)
            print("hudson recombs: ", recombs, " overlapping: ", overlap)

            if print_regions:
                for (a, b) in regions:
                    print("({}, {})".format(a, b))
            continue

        if (not overlap):
            (recombs, regions) = k_hudson_no_overlap(k)

            print(k, "-hudson no overlap recombs: ", recombs)

            if print_regions:
                for (ays, bs) in regions:
                    print("ays: ", ays, "   bs: ", bs)
        elif (k == 2):
            (recombs, regions) = two_hudson()

            print("2-hudson type 1 recombs: ", recombs)

            if print_regions:
                for (a, b, i, j) in regions:
                    print("({}, {}, {}, {})".format(a, b, i, j))

            (recombs, regions) = two_hudson2()

            print("2-hudson type 2 recombs: ", recombs)

            if print_regions:
                for (a, b, i, j) in regions:
                    print("({}, {}, {}, {})".format(a, b, i, j))
        else:
            (recombs, regions) = k_hudson_tail_overlap(k)

            print(k, "-hudson tail overlap recombs: ", recombs)
            print("node that this allows tails to overlap too much")

            if print_regions:
                for (ays, bs) in regions:
                    print("ays: ", ays, "   bs: ", bs)


if __name__ == "__main__":
    main(sys.argv[1:])
