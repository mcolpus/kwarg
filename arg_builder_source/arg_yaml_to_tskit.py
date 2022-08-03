import yaml
import tskit
import pickle
from IPython.display import display, SVG
from simplify_clonal_subtrees import clonal_mrcas

arg = 0
with open("covid_with_clade.yml", "r") as stream:
    try:
        arg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

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
populations = tables.populations

ancestor = 0

# Will treat clades as populations
clades_list = []

for node in arg['nodes']:
    if node['clade'] not in clades_list:
        clades_list.append(node['clade'])

clades_dict = {}
for clade in clades_list:
    id = populations.add_row(pickle.dumps({'name': clade}))
    clades_dict[clade] = id

print(clades_list)
print(clades_dict)




for node in arg["nodes"]:
    if(node["type"] == "root"):
        ancestor = node

    if(node["type"] == "sample"):
        nodes.add_row(flags=1, time=node["time"], population=clades_dict[node['clade']])
    else:
        nodes.add_row(flags=0, time=node["time"], population=clades_dict[node['clade']])

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


muts_id_to_site = {}
for id,site in enumerate(tables.mutations.asdict()["site"]):
    muts_id_to_site[id] = site


tables.build_index()
tables.compute_mutation_times()
ts = tables.tree_sequence()

print(ts.draw_text())


styles = []
# Create a style for each population, programmatically (or just type the string by hand)
for colour, p in zip(['black', 'green', 'blue', 'orange', 'yellow', 'purple', 'red', 'aqua', 'violet', 'cyan', 'gold', 'brown', 'gray'], ts.populations()):
    # target the symbols only (class "sym")
    s = f".node.p{p.id} > .sym " + "{" + f"fill: {colour}" + "}"
    styles.append(s)
    print(f'"{s}" applies to nodes from population {pickle.loads(p.metadata)["name"]} (id {p.id})')
css_string = " ".join(styles)
print(f'CSS string applied:\n    "{css_string}"')


# Height and width for the SVG: optional but useful for this notebook
svg_size = (1600, 500)
svg_string = ts.draw_svg(
    path="svg_vis.svg",
    size=svg_size,
    y_axis=False, y_label=" ",  # optional: show a time scale on the left
    # Match the axis coordinate systems to the text view
    time_scale="rank", x_scale="treewise",
    mutation_labels={}, # muts_id_to_site
    style=css_string,
    node_labels={},
    # x_lim = [0,1000]
)
display(SVG(svg_string))  # If the last line in a cell, wrapping this in display() is not needed


# Now try simplifying

clonal_nodes = clonal_mrcas(ts)
print(f"Subtrees under nodes {clonal_nodes} are clonal")
reduced_ts, node_map = ts.simplify(clonal_nodes, map_nodes=True)

print(
    "The tree seq with clonal subtrees removed",
    "(drawn using the original labels)",
    reduced_ts.draw_text(node_labels = {n:f"{i}" for i, n in enumerate(node_map)}),
    sep="\n",
)

svg_string2 = reduced_ts.draw_svg(
    path="svg_vis2.svg",
    size=svg_size,
    node_labels={n:f"{i}" for i, n in enumerate(node_map)},
    y_axis=False, y_label=" ",  # optional: show a time scale on the left
    # Match the axis coordinate systems to the text view
    time_scale="rank", x_scale="treewise",
    style=css_string,
    mutation_labels={}, # muts_id_to_site
    # x_lim = [0,12000]
)
display(SVG(svg_string2))  # If the last line in a cell, wrapping this in display() is not needed
