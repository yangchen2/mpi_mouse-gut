#!/bin/bash -l
#SBATCH --job-name=tax_assign
#SBATCH --output=/ddn_scratch/yac027/05_gallo_metaG/4_post-qiita/mpi/logs/tax_assign.out
#SBATCH --error=/ddn_scratch/yac027/05_gallo_metaG/4_post-qiita/mpi/logs/tax_assign.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=highmem
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# AUTHOR: Yang Chen
# DATE: 2025-12-09
# DESCRIPTION: Map taxonomy to filtered BIOM tables and create collapsed versions
#              at species and genus levels

source ~/.bashrc
source activate qiime2-2023.2

BASE_DIR="/ddn_scratch/yac027/05_gallo_metaG/4_post-qiita/mpi"
TABLE_PATH="${BASE_DIR}/tables"
TAX_FILE="${BASE_DIR}/reference/WoLr2.tax"

# Create logs directory if it doesn't exist
mkdir -p ${BASE_DIR}/logs

echo "=============================================================="
echo "[$(date)] Starting taxonomy mapping and collapse"
echo "=============================================================="

# Process each 2_ BIOM file
for BIOM_FILE in ${TABLE_PATH}/2_*.biom; do
    if [ ! -f "${BIOM_FILE}" ]; then
        continue
    fi
    
    # Extract base name without extension
    BASENAME=$(basename "${BIOM_FILE}" .biom)
    
    echo "[INFO] Processing ${BASENAME}"
    
    # Export variables for Python
    export BIOM_FILE TAX_FILE TABLE_PATH BASENAME
    
    python <<'PYCODE'
import os
import pandas as pd
import biom
from collections import defaultdict

biom_file = os.environ["BIOM_FILE"]
tax_file = os.environ["TAX_FILE"]
table_path = os.environ["TABLE_PATH"]
basename = os.environ["BASENAME"]

# Load taxonomy file
print(f"[INFO] Loading taxonomy from {tax_file}")
tax_dict = {}
with open(tax_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            genome_id = parts[0]
            taxonomy = parts[1]
            tax_dict[genome_id] = taxonomy

# Load BIOM table
print(f"[INFO] Loading BIOM table: {biom_file}")
table = biom.load_table(biom_file)
df = table.to_dataframe(dense=True)

# Function to extract taxonomic level
def extract_tax_level(tax_string, level):
    """Extract specific taxonomic level from taxonomy string"""
    levels = tax_string.split('; ')
    for lev in levels:
        if lev.startswith(level):
            return lev
    return f"{level}unclassified"

# Map taxonomy and create gOTU versions with unique IDs
print("[INFO] Creating gOTU species and genus tables...")

# Create species and genus mapped DataFrames with unique IDs
species_map = {}
genus_map = {}

for genome_id in df.index:
    if genome_id in tax_dict:
        tax_string = tax_dict[genome_id]
        species = extract_tax_level(tax_string, 's__')
        genus = extract_tax_level(tax_string, 'g__')
        # Keep genome ID to ensure uniqueness
        species_map[genome_id] = f"{genome_id}_gOTU_{species}"
        genus_map[genome_id] = f"{genome_id}_gOTU_{genus}"
    else:
        species_map[genome_id] = f"{genome_id}_gOTU_s__unclassified"
        genus_map[genome_id] = f"{genome_id}_gOTU_g__unclassified"

# Create renamed DataFrames
df_species = df.copy()
df_species.index = df_species.index.map(lambda x: species_map.get(x, f"{x}_gOTU_s__unclassified"))

df_genus = df.copy()
df_genus.index = df_genus.index.map(lambda x: genus_map.get(x, f"{x}_gOTU_g__unclassified"))

# Save gOTU versions
# Species gOTU
species_gotu_file = os.path.join(table_path, f"3_{basename[2:]}_gOTU_species.biom")
species_gotu_table = biom.Table(
    df_species.values,
    observation_ids=df_species.index.astype(str),
    sample_ids=df_species.columns.astype(str)
)
with biom.util.biom_open(species_gotu_file, 'w') as f:
    species_gotu_table.to_hdf5(f, "gOTU species level")
print(f"[INFO] Saved: {species_gotu_file}")

# Genus gOTU
genus_gotu_file = os.path.join(table_path, f"3_{basename[2:]}_gOTU_genus.biom")
genus_gotu_table = biom.Table(
    df_genus.values,
    observation_ids=df_genus.index.astype(str),
    sample_ids=df_genus.columns.astype(str)
)
with biom.util.biom_open(genus_gotu_file, 'w') as f:
    genus_gotu_table.to_hdf5(f, "gOTU genus level")
print(f"[INFO] Saved: {genus_gotu_file}")

# Create collapsed versions
print("[INFO] Creating collapsed species and genus tables...")

# Collapse to species
species_collapsed = defaultdict(lambda: pd.Series(0, index=df.columns))
for genome_id in df.index:
    if genome_id in tax_dict:
        species = extract_tax_level(tax_dict[genome_id], 's__')
    else:
        species = "s__unclassified"
    species_collapsed[species] += df.loc[genome_id]

df_species_collapsed = pd.DataFrame(species_collapsed).T

# Collapse to genus
genus_collapsed = defaultdict(lambda: pd.Series(0, index=df.columns))
for genome_id in df.index:
    if genome_id in tax_dict:
        genus = extract_tax_level(tax_dict[genome_id], 'g__')
    else:
        genus = "g__unclassified"
    genus_collapsed[genus] += df.loc[genome_id]

df_genus_collapsed = pd.DataFrame(genus_collapsed).T

# Save collapsed versions
# Species collapsed
species_collapsed_file = os.path.join(table_path, f"3_{basename[2:]}_species_collapsed.biom")
species_collapsed_table = biom.Table(
    df_species_collapsed.values,
    observation_ids=df_species_collapsed.index.astype(str),
    sample_ids=df_species_collapsed.columns.astype(str)
)
with biom.util.biom_open(species_collapsed_file, 'w') as f:
    species_collapsed_table.to_hdf5(f, "species collapsed")
print(f"[INFO] Saved: {species_collapsed_file}")

# Genus collapsed
genus_collapsed_file = os.path.join(table_path, f"3_{basename[2:]}_genus_collapsed.biom")
genus_collapsed_table = biom.Table(
    df_genus_collapsed.values,
    observation_ids=df_genus_collapsed.index.astype(str),
    sample_ids=df_genus_collapsed.columns.astype(str)
)
with biom.util.biom_open(genus_collapsed_file, 'w') as f:
    genus_collapsed_table.to_hdf5(f, "genus collapsed")
print(f"[INFO] Saved: {genus_collapsed_file}")

print(f"[INFO] Completed processing {basename}")
print("=" * 60)
PYCODE

done

echo "[$(date)] All taxonomy mapping and collapse completed"
echo "=============================================================="
