import networkx as nx
import matplotlib.pylab as plt
import sys

if len(sys.argv) != 2:
    print("Should input filename only")
else:
    filename = sys.argv[1]
    print(filename)
    g = nx.read_gml(filename)
    nx.draw(g)
    plt.savefig(filename[:-4] + ".png")