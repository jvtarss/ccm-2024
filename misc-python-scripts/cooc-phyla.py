import pandas as pd

feature_table = pd.read_csv("gg-phyla-feature-table-hellinger.tsv", sep="\t", index_col=0)
metadata = pd.read_csv("metadata_backup.tsv", sep="\t")

common_samples = set(feature_table.columns).intersection(metadata['sample-id'])
metadata = metadata[metadata['sample-id'].isin(common_samples)]

condition_groups = metadata['condition'].unique()

phyla_by_group = {group: set() for group in condition_groups}

for group in condition_groups:
    samples_in_group = metadata[metadata['condition'] == group]['sample-id']
    
    group_feature_table = feature_table[samples_in_group]
    
    phyla_in_group = group_feature_table[group_feature_table > 0].dropna(how='all').index
    phyla_in_group = [phylum.split(';')[-1].replace('p__', '') for phylum in phyla_in_group]
    
    phyla_by_group[group].update(phyla_in_group)

common_phyla = set.intersection(*phyla_by_group.values())

with open("cooc-phyla.txt", "w") as f:
    for phylum in sorted(common_phyla):
        f.write(f"{phylum}\n")

print(f"arquivo 'cooc-phyla.txt' criado com {len(common_phyla)} filos comuns.")
