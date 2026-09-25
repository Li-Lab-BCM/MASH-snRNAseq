# =============================================================================
# 02_merge.R
# Merge the 13 QC-filtered samples and attach sample metadata.
#
# Input : <DATA_DIR>/<Sample>-cellbender-SCDBL_final.rds
# Output: <DATA_DIR>/2026MASH_FINAL_merged_QC1.rds
# =============================================================================

library(Seurat)

DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")

meta <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))

objs <- lapply(seq_len(nrow(meta)), function(i) {
  obj <- readRDS(file.path(DATA_DIR, paste0(meta$Sample[i], "-cellbender-SCDBL_final.rds")))
  obj$sample     <- meta$Sample[i]
  obj$orig.ident <- meta$orig.ident[i]
  obj$group      <- meta$Group[i]
  obj$replicate  <- meta$Replicate[i]
  obj
})

seurat_object <- merge(objs[[1]], y = objs[-1],
                       add.cell.ids = meta$orig.ident,
                       project = "MergedSamples")
rm(objs)

message("Total cells: ", ncol(seurat_object))
print(table(seurat_object$orig.ident))

saveRDS(seurat_object, file.path(DATA_DIR, "2026MASH_FINAL_merged_QC1.rds"))

sessionInfo()
