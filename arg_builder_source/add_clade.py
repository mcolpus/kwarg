import yaml
import pandas
 
# reading the CSV file
csvFile = pandas.read_csv('metadata.csv')
 
# strain_to_clade = csvFile.loc[:, ["strain", "GISAID_clade"]]

strain_to_clade = {}
for idx, row in csvFile.iterrows():
    strain_to_clade[row.strain] = row.GISAID_clade


arg = 0
with open("covid.yml", "r") as stream:
    try:
        arg = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

for node in arg['nodes']:
    label = node['label']
    clade = '?'
    if 'hCoV-19' in label:
        clade = strain_to_clade[label]
    node['clade'] = clade


with open('covid_with_clade.yml', 'w') as file:
    outputs = yaml.dump(arg, file)
    