#!/bin/bash

#SBATCH -n 10
#SBATCH --partition="YOUR_PARTITION"
#SBATCH --mem-per-cpu=30000M
#SBATCH --output=output.out
#SBATCH --error=error.err

# ==============================
# USER INPUTS (EDIT THESE ONLY)
# ==============================

# Path to scDRS script
SCDRS_SCRIPT="/path/to/scDRS/compute_downstream.py"

# Input data
H5AD_FILE="/path/to/input_file.h5ad"
SCORE_FILE="/path/to/score_file.gz"

# Output directory
OUT_DIR="/path/to/output_directory"

# Analysis parameters
CELL_TYPE="development_stage"

# ==============================
# RUN
# ==============================

python "$SCDRS_SCRIPT" \
  --h5ad_file "$H5AD_FILE" \
  --score_file "$SCORE_FILE" \
  --cell_type "$CELL_TYPE" \
  --flag_gene True \
  --flag_filter False \
  --flag_raw_count True \
  --out_folder "$OUT_DIR"