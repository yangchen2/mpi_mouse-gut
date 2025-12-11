#!/bin/bash -l
#SBATCH --job-name=micov_filt
#SBATCH --output=/ddn_scratch/yac027/05_gallo_metaG/4_post-qiita/mpi/logs/micov_filt.out
#SBATCH --error=/ddn_scratch/yac027/05_gallo_metaG/4_post-qiita/mpi/logs/micov_filt.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=highmem
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=yac027@ucsd.edu

# AUTHOR: Yang Chen
# DATE: 2025-12-09
# DESCRIPTION: Filter BIOM tables based on micov percent_covered > threshold (%).
#              Running thresholds: 1%,  5%, and 10%.

source ~/.bashrc
source activate qiime2-2023.2

BASE_DIR="/ddn_scratch/yac027/05_gallo_metaG/4_post-qiita/mpi"
MICOV_FILE="${BASE_DIR}/coverage/coverage_calculation.tsv"
BIOM_IN="${BASE_DIR}/tables/1_228413_gOTU.biom"
TABLE_PATH="${BASE_DIR}/tables"

# Thresholds to evaluate (percent)
THRESHOLDS=(1 5 10)

# Create logs directory if it doesn't exist
mkdir -p ${BASE_DIR}/logs

echo "=============================================================="
echo "[$(date)] Starting micov_filt"
echo "=============================================================="

# Sanity checks
if [[ ! -f "${MICOV_FILE}" ]]; then
    echo "[ERROR] Missing micov file: ${MICOV_FILE}"
    exit 1
fi

if [[ ! -f "${BIOM_IN}" ]]; then
    echo "[ERROR] Missing BIOM file: ${BIOM_IN}"
    exit 1
fi

# Export inputs to Python
export BIOM_IN MICOV_FILE TABLE_PATH

# Run Python once and handle all thresholds inside
python <<'PYCODE'
import os
import pandas as pd
import biom

biom_in = os.environ["BIOM_IN"]
micov_file = os.environ["MICOV_FILE"]
table_path = os.environ["TABLE_PATH"]

# Thresholds from shell
THRESHOLDS = [1, 5, 10]

print(f"[INFO] Loading BIOM table: {biom_in}")
table = biom.load_table(biom_in)
df = table.to_dataframe(dense=True)

# Load micov coverage table
micov = pd.read_csv(micov_file, sep="\t")
micov['genome_id'] = micov['genome_id'].astype(str)

available_ids = set(df.index.astype(str))

for thresh in THRESHOLDS:
    print(f"\n[INFO] Applying micov threshold: >{thresh}%")

    keep_ids = set(
        micov.loc[micov['percent_covered'] > float(thresh), 'genome_id'].astype(str)
    )

    matched = keep_ids.intersection(available_ids)
    filtered_df = df.loc[df.index.isin(matched)]

    print(f"[INFO] {len(matched)} genome_ids matched out of {len(keep_ids)} above {thresh}%.")
    print(f"[INFO] Filtered table size: {filtered_df.shape[0]} features x {filtered_df.shape[1]} samples.")

    biom_out = os.path.join(table_path, f"2_228413_gOTU_micov{thresh}%.biom")

    # Convert back to BIOM
    filtered_table = biom.Table(
        filtered_df.values,
        observation_ids=filtered_df.index.astype(str),
        sample_ids=filtered_df.columns.astype(str)
    )

    with biom.util.biom_open(biom_out, 'w') as f:
        filtered_table.to_hdf5(f, f"filtered biom >{thresh}%")

    print(f"[INFO] Saved filtered BIOM: {biom_out}")

print("\n[INFO] All thresholds complete.")
PYCODE

echo "[$(date)] Finished micov_filt"
echo "=============================================================="
