import yaml
import msprime
import tskit
from IPython.display import display, SVG
import numpy as np

def clonal_mrcas(ts):
    # full length nodes without recombination will not be a parent in the edges in or out, apart from at the start
    for interval, edges_out, edges_in in ts.edge_diffs():
        if interval.left==0:
            is_full_length_clonal = np.ones(ts.num_nodes, dtype=bool)  # all nodes start as clonal
        else:
            for e in edges_in:
                is_full_length_clonal[e.parent] = False
        for e in edges_out:
            is_full_length_clonal[e.parent] = False
    clonal_nodes = np.where(is_full_length_clonal)[0]

    tables = ts.dump_tables()
    edges = tables.edges
    # only keep edges where both the child and the parent are full-length clonal nodes
    keep_edge = np.logical_and(np.isin(edges.child, clonal_nodes),  np.isin(edges.parent, clonal_nodes))
    tables.edges.set_columns(
            left = tables.edges.left[keep_edge],
            right=tables.edges.right[keep_edge],
            parent=tables.edges.parent[keep_edge],
            child=tables.edges.child[keep_edge],
        )
    ts = tables.tree_sequence()
    print("Debug: show the clonal subtrees", ts.draw_text(), sep="\n")
    assert ts.num_trees == 1
    return ts.first().roots

ts = msprime.simulate(8, recombination_rate=1, random_seed=321)
print("Original ts", ts.draw_text(), sep="\n")

clonal_nodes = clonal_mrcas(ts)
print(f"Subtrees under nodes {clonal_nodes} are clonal")
reduced_ts, node_map = ts.simplify(clonal_nodes, map_nodes=True)

print(
    "The tree seq with clonal subtrees removed",
    "(drawn using the original labels)",
    reduced_ts.draw_text(node_labels = {n:f"{i}" for i, n in enumerate(node_map)}),
    sep="\n",
)