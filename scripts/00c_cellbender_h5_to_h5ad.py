"""
00c_cellbender_h5_to_h5ad.py
Convert CellBender filtered output (.h5) to .h5ad for reading into R
with zellkonverter (scripts/01_QC_doublets.R).

Input : <DATA_DIR>/<Sample>/cellbender_cleaned_filtered.h5
Output: <DATA_DIR>/<Sample>/cellbender_cleaned_filtered.h5ad
"""
import csv
import os

import numpy as np
import scanpy as sc
import scipy.sparse

DATA_DIR = "/path/to/MASHproject"
REPO_DIR = os.path.join(DATA_DIR, "MASH-snRNAseq")

remove_zero_genes = True  # drop genes with zero counts in all cells

with open(os.path.join(REPO_DIR, "metadata", "sample_metadata.tsv")) as f:
    samples = [r["Sample"] for r in csv.DictReader(f, delimiter="\t")]

for sample in samples:
    input_path = os.path.join(DATA_DIR, sample, "cellbender_cleaned_filtered.h5")
    output_path = os.path.join(DATA_DIR, sample, "cellbender_cleaned_filtered.h5ad")

    adata = sc.read_10x_h5(input_path)
    adata.var_names_make_unique()

    if remove_zero_genes:
        nonzero = np.array((adata.X > 0).sum(axis=0)).flatten() > 0
        adata = adata[:, nonzero]

    if not scipy.sparse.issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)

    adata.raw = adata
    adata.write(output_path)
    print(f"{sample}: saved {output_path}, shape: {adata.shape}")
