import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('phylum-table.tsv', sep='\t', skiprows=2, index_col=0)

def hellinger_transform(df):
    relative_abundance = df.div(df.sum(axis=0), axis=1)
    hellinger_df = np.sqrt(relative_abundance)
    return hellinger_df

hellinger_data = hellinger_transform(data)

hellinger_data.to_csv('hellinger_transformed_data.csv')

abundance_sum = hellinger_data.sum(axis=1)

most_abundant_phyla = abundance_sum.sort_values(ascending=False)

most_abundant_phyla.to_csv('most_abundant_phyla.csv')

