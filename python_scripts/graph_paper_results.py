import numpy as np
import sys
import getopt
import matplotlib as mpl
from matplotlib import pyplot as plt
from graph_results import read_args, calculate_optimal_curve

np.set_printoptions(edgeitems=30, linewidth=100000,
                    formatter=dict(float=lambda x: "%.3g" % x))

plt.rc('font', size=12) 

# Reads from csv containing many runs and create a grid for optimal args


def graph_all(inputfile, num_files, outputfile):
    fig = plt.figure(1, figsize=[8, 6])

    for i in range(num_files+1):
        file = inputfile + "_" + str(i) + ".csv"
        print(file)
        arg_runs = read_args(file)
        (fewest_rare_muts, min_recombs, max_recombs, recombinations, rare_muts) = calculate_optimal_curve(arg_runs, True)

        print(fewest_rare_muts)


        plt.plot(range(min_recombs, max_recombs+1), fewest_rare_muts[min_recombs:],
                color='red', alpha=0.1, label="optimal")



    plt.title('Recombinations vs Recurrent mutations \n for all samples')
    plt.xlabel('Recombinations')
    plt.ylabel('recurrent mutations (all kinds)')

    # plt.show()

    fig.savefig(outputfile, dpi=200, bbox_inches='tight')


def calculate_statistics(inputfile, num_files, outputfile):
    samples_stats = list()

    for i in range(num_files+1):
        file = inputfile + "_" + str(i) + ".csv"
        arg_runs = read_args(file)
        (fewest_rare_muts, min_recombs, max_recombs, recombinations, rare_muts) = calculate_optimal_curve(arg_runs, True)

        only_rms_min = fewest_rare_muts[0]
        only_recomb_min = max_recombs
        both_min = 1000000
        for j in range(min_recombs, max_recombs+1):
            both_min = min(both_min, j + fewest_rare_muts[j])
        
        samples_stats.append((only_rms_min, only_recomb_min, both_min))


    only_rms = np.array([rm for rm, recomb, both in samples_stats])
    q1_rm, q2_rm, q3_rm = np.percentile(only_rms, [25, 50, 75], interpolation='midpoint')
    print(q1_rm, q2_rm, q3_rm)
    only_recombs = np.array([recomb for rm, recomb, both in samples_stats])
    q1_recomb, q2_recomb, q3_recomb = np.percentile(only_recombs, [25, 50, 75], interpolation='midpoint')
    print(q1_recomb, q2_recomb, q3_recomb)
    boths = np.array([both for rm, recomb, both in samples_stats])
    q1_both, q2_both, q3_both = np.percentile(boths, [25, 50, 75], interpolation='midpoint')
    print(q1_both, q2_both, q3_both)

    fig, ax = plt.subplots()
    ax.set_title('Hide Outlier Points')
    ax.boxplot([only_rms, only_recombs, boths], showfliers=True, labels=["min rms", "min recombs", "min combined"])

    fig.savefig(outputfile, dpi=200, bbox_inches='tight')

    return samples_stats

def calculate_all_statistics(inputfile, num_files, outputfile):
    stats_by_group = list()
    num_groups = len(inputfile)

    for input, num in zip(inputfile, num_files):

        rms = list()
        recombs = list()
        boths = list()

        for i in range(num+1):
            file = input + "_" + str(i) + ".csv"
            arg_runs = read_args(file)
            (fewest_rare_muts, min_recombs, max_recombs, _, _) = calculate_optimal_curve(arg_runs, True)

            rms.append(fewest_rare_muts[0])
            recombs.append(max_recombs)

            both_min = 1000000
            for j in range(min_recombs, max_recombs+1):
                both_min = min(both_min, j + fewest_rare_muts[j])

            boths.append(both_min)
        
        stats_by_group.append((rms, recombs, boths))

    all_stats = [l for stats in stats_by_group for l in stats]
    positions = [i*0.6 + (i//3) for i in range(len(all_stats))]
    colors = ["red", "blue", "green"] * num_groups
    labels = ['', '', ''] * num_groups
    for i in range(num_groups):
        labels[3*i + 1] = str((i+1)*20)

    fig, ax = plt.subplots()
    ax.set_title('Rare events by sample size')
    bplot = ax.boxplot(all_stats, showfliers=True, positions=positions, patch_artist=True, labels=labels)

    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    fig.savefig(outputfile, dpi=200, bbox_inches='tight')


def main(argv):
    inputfile = ''
    num_files = 0
    outputfile = ''
    calculate_stats = False


    try:
        opts, args = getopt.getopt(argv, "hi:o:n:c", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('graph_results.py -i <input file> -o <output file> -n <number of input files> -c (calculate stats)')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('graph_results.py -i <input file> -o <output file> -n <number of input files> -c (calculate stats)')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-n"):
            num_files = arg
        elif opt in ("-c"):
            calculate_stats = True

    print('Input file is: ', inputfile)
    print('Number of files: ', num_files)
    print('Output file is: ', outputfile)


    if calculate_stats:
        if ',' in inputfile:
            inputfile = inputfile.split(',')
            num_files = [int(n) for n in num_files.split(',')]
            calculate_all_statistics(inputfile, num_files, outputfile)
        else:
            calculate_statistics(inputfile, num_files, outputfile)
    else:
        graph_all(inputfile, num_files, outputfile)




if __name__ == "__main__":
    main(sys.argv[1:])
