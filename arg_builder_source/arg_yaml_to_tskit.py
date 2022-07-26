import yaml
import tskit
from IPython.display import display, SVG

arg = 0
with open("testyaml.yml", "r") as stream:
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

n = arg["number of nodes"]
sequence_length = arg["sequence length"]
print("length", sequence_length)

# individuals = tskit.IndividualTable()
tables = tskit.TableCollection()
tables.sequence_length = sequence_length
nodes = tables.nodes
edges = tables.edges
sites = tables.sites
mutations = tables.mutations

ancestor = 0


for node in arg["nodes"]:
    if(node["type"] == "root"):
        ancestor = node

    if(node["type"] == "sample"):
        nodes.add_row(flags=1, time=node["time"])
    else:
        nodes.add_row(flags=0, time=node["time"])

if(ancestor == 0):
    print("Error, no ancestor node given")
    quit()

ancestral_muts = [int(i) for i in str(ancestor["mutations"]).split(',')] if ancestor["mutations"] != None else []
for i in range(sequence_length):
    sites.add_row(i, '1' if i in ancestral_muts else '0')


def add_edge(start, end, source, target):
    edges.add_row(start, end, source, target)

for edge in arg["edges"]:
    if(edge["type"] == "single"):
        add_edge(0, sequence_length, edge["from"], edge["to"])
    elif(edge["type"] == "prefix"):
        add_edge(0, edge["crossover"], edge["from"], edge["to"])
    elif(edge["type"] == "suffix"):
        add_edge(edge["crossover"], sequence_length, edge["from"], edge["to"])

    if edge["mutations"] != None:
        muts = [int(i) for i in str(edge["mutations"]).split(',')]

        for mut in muts:
            derived_state = 0 if mut in ancestral_muts else 1
            mutations.add_row(mut, edge["to"], str(derived_state))
    
    if edge["back_mutations"] != None:
        muts = [int(i) for i in str(edge["back_mutations"]).split(',')]

        for mut in muts:
            derived_state = 1 if mut in ancestral_muts else 0
            mutations.add_row(mut, edge["to"], str(derived_state))


tables.sort()

print(tables)

muts_id_to_site = {}
for id,site in enumerate(tables.mutations.asdict()["site"]):
    muts_id_to_site[id] = site

print(muts_id_to_site)

tables.build_index()
tables.compute_mutation_times()
ts = tables.tree_sequence()

print(ts.draw_text())

# Height and width for the SVG: optional but useful for this notebook
svg_size = (1600, 500)
svg_string = ts.draw_svg(
    path="svg_vis.svg",
    size=svg_size,
    y_axis=False, y_label=" ",  # optional: show a time scale on the left
    # Match the axis coordinate systems to the text view
    time_scale="rank", x_scale="treewise",
    mutation_labels=muts_id_to_site
)
display(SVG(svg_string))  # If the last line in a cell, wrapping this in display() is not needed
