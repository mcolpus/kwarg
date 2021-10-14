import yaml
import tskit
from IPython.display import display, SVG

arg = 0
with open("example_ts.yml", "r") as stream:
    try:
        arg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

# header: arg exported to text
# number_of_nodes: 27
# nodes:
#   - id: 0
#     label: 1
#     sequence: 0010000001000001001101110111101010101000000
#     type: sample
# edges:
#   - source: 14
#     target: 0
#     type: single
#     position: 0
#     mutations:
#       - 9
#       - -8

n = arg["number_of_nodes"]
sequence_length = arg["sequence_length"]
print("length", sequence_length)

# individuals = tskit.IndividualTable()
tables = tskit.TableCollection()
tables.sequence_length = sequence_length
nodes = tables.nodes
edges = tables.edges
sites = tables.sites
mutations = tables.mutations

ancestor = 0

# Note that the nodes are given in reverse order due to the way that coalescence algorithm works
time = 1
for node in arg["nodes"]:
    if(node["type"] == "ancestor"):
        ancestor = node

    if(node["type"]=="sample"):
        nodes.add_row(flags=1, time = 0)
    else:
        nodes.add_row(flags=0,time=time)
        time += 1

if(ancestor == 0):
    print("Error, no ancestor node given")
    quit()
for i in range(sequence_length):
    sites.add_row(i, ancestor["sequence"][i])

def add_edge(start, end, source, target):
    edges.add_row(start, end, source, target)

for edge in arg["edges"]:
    if(edge["type"] == "single"):
        add_edge(0, sequence_length, edge["source"], edge["target"])
        # edges.add_row(0, sequence_length, edge["source"], edge["target"])
    elif(edge["type"] == "prefix"):
        add_edge(0, edge["position"], edge["source"], edge["target"])
    elif(edge["type"] == "suffix"):
        add_edge(edge["position"], sequence_length, edge["source"], edge["target"])

    if "mutations" in edge:
        for mutation in edge["mutations"]:
            derived_state = arg["nodes"][edge["target"]]["sequence"][mutation]
            mutations.add_row(mutation, edge["target"], derived_state)


print(tables)

tables.sort()
tables.build_index()
tables.compute_mutation_times()
ts = tables.tree_sequence()

print(ts.draw_text())

svg_size = (1600, 500) # Height and width for the SVG: optional but useful for this notebook
svg_string = ts.draw_svg(
    path="svg_vis.svg",
    size=svg_size,
    y_axis=True, y_label=" ",  # optional: show a time scale on the left
    time_scale="rank", x_scale="treewise",  # Match the axis coordinate systems to the text view
)
# display(SVG(svg_string))  # If the last line in a cell, wrapping this in display() is not needed