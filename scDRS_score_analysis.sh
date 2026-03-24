#!/bin/bash

#SBATCH -n 10
#SBATCH --partition="YOUR_PARTITION"
#SBATCH --mem-per-cpu=100000M
#SBATCH --output=output.out
#SBATCH --error=error.err

# -----------------------------
# User-defined paths and inputs
# -----------------------------
SCDRS_SCRIPT="/path/to/scDRS/compute_score.py"
H5AD_FILE="/path/to/input_file.h5ad"
GS_FILE="/path/to/Gene_set.file"
OUT_FOLDER="/path/to/output_folder"

# -----------------------------
# Fixed parameters
# -----------------------------
H5AD_SPECIES="human"
GS_SPECIES="human"
N_CTRL=1000

python "$SCDRS_SCRIPT" \
  --h5ad_file "$H5AD_FILE" \
  --h5ad_species "$H5AD_SPECIES" \
  --gs_file "$GS_FILE" \
  --gs_species "$GS_SPECIES" \
  --flag_filter True \
  --flag_raw_count True \
  --flag_return_ctrl_raw_score True \
  --flag_return_ctrl_norm_score True \
  --out_folder "$OUT_FOLDER" \
  --n_ctrl "$N_CTRL"