#!/usr/bin/env bash
# CellBender remove-background for all 13 snRNA-seq libraries.
# CellBender version: 0.3.0
# Unlisted parameters use CellBender defaults (e.g. --fpr 0.01).

set -euo pipefail

DATA_DIR="${MASH_DATA_DIR:-/path/to/MASHproject}"

run_cellbender() {
  sample="$1"
  expected_cells="$2"

  sample_dir="${DATA_DIR}/${sample}"
  input_file="${sample_dir}/raw_feature_bc_matrix.h5"
  output_file="${sample_dir}/cellbender_cleaned.h5"

  if [[ ! -f "${input_file}" ]]; then
    echo "Input file not found: ${input_file}" >&2
    exit 1
  fi

  echo "Running CellBender for ${sample}"

  cellbender remove-background \
    --input "${input_file}" \
    --output "${output_file}" \
    --expected-cells "${expected_cells}" \
    --total-droplets-included 40000 \
    --cuda \
    --z-dim 50 \
    --epochs 100 \
    --learning-rate 1e-4 \
    --posterior-batch-size 64 \
    --low-count-threshold 10
}

# CHOW
run_cellbender A1 20000
run_cellbender A2 18000
run_cellbender A3 16000

# MASH
run_cellbender B1 20000
run_cellbender B2 22000
run_cellbender B3 14500

# MASH + anti-OPN
run_cellbender D1 16000
run_cellbender D2 12300
run_cellbender D3 14500

# MASH Spp1 knockout
run_cellbender F1 15000
run_cellbender F3 13600

# MASH flox
run_cellbender G1 11500
run_cellbender G2 16500

echo "CellBender processing completed."
