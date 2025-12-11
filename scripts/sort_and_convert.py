import pandas as pd
import biom
from biom import load_table

# Load in biom and metadata
taxon_level = 'genus'

biom_path = f"../tables/3_228413_gOTU_micov75%_{taxon_level}_collapsed.biom"
metadata = pd.read_csv("../metadata/marta_metadata.tsv", sep='\t')

# Organize index names
biom_table = load_table(biom_path)
df = biom_table.to_dataframe()  # BIOM already has features as rows, samples as columns
df.columns = df.columns.astype(str)
df.columns = df.columns.str.replace('16207.', '')
df = df[metadata['sample_name']]
df = df.transpose()

# Save this as tsv table
tsv_table_path = f"../tables/3_228413_gOTU_micov75%_{taxon_level}_collapsed.tsv"
df.to_csv(tsv_table_path, sep='\t')
