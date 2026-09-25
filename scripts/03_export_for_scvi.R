# =============================================================================
# 03_export_for_scvi.R
# Export the merged QC1 object to .h5ad (raw counts) for scVI in Python.
#
# Input : <DATA_DIR>/2026MASH_FINAL_merged_QC1.rds
# Output: <DATA_DIR>/seurat_for_scvi_QC1.h5ad
# =============================================================================


library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)

DATA_DIR <- "/path/to/MASHproject"

obj <- readRDS(file.path(DATA_DIR, "2026MASH_FINAL_merged_QC1.rds"))
obj <- JoinLayers(obj)

sce <- as.SingleCellExperiment(obj, assay = "RNA")
writeH5AD(sce, file.path(DATA_DIR, "seurat_for_scvi_QC1.h5ad"), X_name = "counts")

message("Exported ", ncol(sce), " cells x ", nrow(sce), " genes")
sessionInfo()
