import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv('phylum-table.tsv', sep='\t', skiprows=2, index_col=0)

# Perform Hellinger transformation
def hellinger_transform(df):
    # Calculate relative abundance
    relative_abundance = df.div(df.sum(axis=0), axis=1)
    # Apply square root transformation
    hellinger_df = np.sqrt(relative_abundance)
    return hellinger_df

hellinger_data = hellinger_transform(data)

# Save the transformed data
hellinger_data.to_csv('hellinger_transformed_data.csv')

# Sum the abundances across all samples
abundance_sum = hellinger_data.sum(axis=1)

# Sort phyla by abundance
most_abundant_phyla = abundance_sum.sort_values(ascending=False)

# Save the most abundant phyla
most_abundant_phyla.to_csv('most_abundant_phyla.csv')

# Plot the top 20 most abundant phyla
top_20_phyla = most_abundant_phyla.head(20)
plt.figure(figsize=(10, 8))
top_20_phyla.plot(kind='bar', color='skyblue')
plt.title('Top 20 Most Abundant Phyla')
plt.xlabel('Phylum')
plt.ylabel('Abundance (Hellinger Transformed)')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('top_20_abundant_phyla.png')
plt.show()
