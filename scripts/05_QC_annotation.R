# =============================================================================
# 05_QC_annotation.R
# Clustering on the scVI latent space, cell-type annotation and removal of
# contaminating / low-quality cells. Final object uses Harmony for
# visualisation and downstream analyses.
#
# Input : <DATA_DIR>/2026MASH_FINAL_merged_QC1.rds
#         <DATA_DIR>/seurat_scvi_5group.h5ad        (04_scvi_integration.py)
# Output: <DATA_DIR>/seurat_object_QC4.rds          (final object)
#         logs/cleanup/cell_filtering_steps.csv
#         metadata/cell_annotations.csv.gz
# =============================================================================

library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)
library(harmony)
library(dplyr)

DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")
dir.create(file.path(REPO_DIR, "logs", "cleanup"), recursive = TRUE, showWarnings = FALSE)
set.seed(1234)

steps <- data.frame()
log_step <- function(obj, step) {
  steps <<- rbind(steps, data.frame(step = step, n_cells = ncol(obj)))
  message(sprintf("%-40s %d cells", step, ncol(obj)))
}
drop_cells <- function(obj, cells) subset(obj, cells = setdiff(colnames(obj), cells))

# -----------------------------------------------------------------------------
# 1. Load merged object and add scVI latent space
# -----------------------------------------------------------------------------
seurat_object <- readRDS(file.path(DATA_DIR, "2026MASH_FINAL_merged_QC1.rds"))
seurat_object <- JoinLayers(seurat_object)

sce <- readH5AD(file.path(DATA_DIR, "seurat_scvi_5group.h5ad"))
emb <- as.matrix(reducedDim(sce, "X_scVI"))
rownames(emb) <- colnames(sce)
emb <- emb[colnames(seurat_object), ]
colnames(emb) <- paste0("scvi_", seq_len(ncol(emb)))
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = emb, key = "scvi_", assay = "RNA")
rm(sce, emb)
log_step(seurat_object, "Merged")

# -----------------------------------------------------------------------------
# 2. Hemoglobin filter (percent.rb = % Hba/Hbb reads)
# -----------------------------------------------------------------------------
seurat_object <- subset(seurat_object, subset = percent.rb <= 0.2)
log_step(seurat_object, "Hemoglobin filter")

# -----------------------------------------------------------------------------
# 3. Clustering (scVI) and annotation; remove contaminating / low-quality
#    clusters identified from marker genes
# -----------------------------------------------------------------------------
seurat_object <- RunUMAP(seurat_object, reduction = "scvi", dims = 1:30,
                         reduction.name = "umap_scvi") %>%
  FindNeighbors(reduction = "scvi", dims = 1:30) %>%
  FindClusters(resolution = 0.5,algorithm = 1)

seurat_object <- subset(seurat_object, seurat_clusters %in% c(8, 22, 27, 28), invert = TRUE)

new_names <- c(
  "0"  = "Endocardial cell",   "1"  = "Endocardial cell",
  "2"  = "Epicardial cell",    "3"  = "Fibroblast",
  "4"  = "Fibroblast",         "5"  = "Endothelial cell",
  "6"  = "MP/DC",              "7"  = "Atrial Cardiomyocyte",
  "9"  = "MP/DC",              "10" = "Mural cell",
  "11" = "MP/DC",              "12" = "Endocardial cell",
  "13" = "T/NK cell",          "14" = "MP/DC",
  "15" = "LEC",                "16" = "Endocardial cell",
  "17" = "B cell",             "18" = "Proliferating cell",
  "19" = "Epicardial cell",    "20" = "T/NK cell",
  "21" = "Endothelial cell",   "23" = "Schwann cell",
  "24" = "T/NK cell",          "25" = "Adipocyte",
  "26" = "MP/DC"
)
Idents(seurat_object) <- "seurat_clusters"
seurat_object <- RenameIdents(seurat_object, new_names)
seurat_object$celltype <- as.character(Idents(seurat_object))

# Sub-cluster selected cell types and remove small sub-clusters co-expressing
# markers of other lineages and separated from the main population on UMAP
subcluster_ids <- function(obj, ct) {
  sub <- subset(obj, subset = celltype == ct) %>%
    FindNeighbors(reduction = "scvi", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.5,algorithm = 1,verbose = FALSE)
  setNames(as.character(sub$seurat_clusters), colnames(sub))
}
exclude <- list("Epicardial cell" = c("5", "6"),
                "Fibroblast"      = "6",
                "MP/DC"           = c("4", "6", "7"))
cells_remove <- unlist(lapply(names(exclude), function(ct) {
  ids <- subcluster_ids(seurat_object, ct)
  names(ids)[ids %in% exclude[[ct]]]
}))
seurat_object <- drop_cells(seurat_object, cells_remove)
log_step(seurat_object, "Contaminating clusters removed")

# -----------------------------------------------------------------------------
# 4. Harmony integration
# -----------------------------------------------------------------------------
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunHarmony("orig.ident")
seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, min.dist = 0.5)

# -----------------------------------------------------------------------------
# 5. Proliferating cells: assign to nearest cell-type centroid (scVI UMAP);
#    keep those assigned to MP/DC or endothelial cells
# -----------------------------------------------------------------------------
umap_coords <- Embeddings(seurat_object, "umap_scvi")
ct <- seurat_object$celltype
is_prolif <- ct == "Proliferating cell"
centers <- aggregate(umap_coords[!is_prolif, ], by = list(celltype = ct[!is_prolif]), FUN = mean)
rownames(centers) <- centers$celltype
centers <- as.matrix(centers[, -1])
for (i in which(is_prolif)) {
  ct[i] <- names(which.min(sqrt(colSums((t(centers) - umap_coords[i, ])^2))))
}
remove_prolif <- colnames(seurat_object)[is_prolif & !ct %in% c("MP/DC", "Endothelial cell")]
seurat_object$celltype <- ct
seurat_object <- drop_cells(seurat_object, remove_prolif)
log_step(seurat_object, "Proliferating cells reassigned")

# -----------------------------------------------------------------------------
# 6. Remove non-cardiomyocytes with cardiomyocyte transcript signal
# -----------------------------------------------------------------------------
seurat_object <- AddModuleScore(seurat_object, name = "CM_contam_score",
                                features = list(c("Myl7", "Nppa", "Myh6", "Tnnt2", "Tnnc1", "Actc1")))
seurat_object <- drop_cells(seurat_object, colnames(seurat_object)[
  seurat_object$celltype != "Atrial Cardiomyocyte" & seurat_object$CM_contam_score1 > 1])
log_step(seurat_object, "CM-contaminated cells removed")

# -----------------------------------------------------------------------------
# 7. Remove UMAP outliers (> 3 MAD from the cell-type median)
# -----------------------------------------------------------------------------
umap_coords <- Embeddings(seurat_object, "umap")
is_outlier <- rep(FALSE, ncol(seurat_object))
for (x in unique(seurat_object$celltype)) {
  idx <- which(seurat_object$celltype == x)
  is_outlier[idx] <-
    abs(umap_coords[idx, 1] - median(umap_coords[idx, 1])) > 3 * mad(umap_coords[idx, 1]) |
    abs(umap_coords[idx, 2] - median(umap_coords[idx, 2])) > 3 * mad(umap_coords[idx, 2])
}
seurat_object <- seurat_object[, !is_outlier]
log_step(seurat_object, "UMAP outliers removed")

seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, min.dist = 0.5)

# -----------------------------------------------------------------------------
# 8. Final removal of contaminating MP/DC and fibroblast sub-clusters
# -----------------------------------------------------------------------------
recluster <- function(obj, res) {
  NormalizeData(obj) %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 20, verbose = FALSE) %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(reduction = "harmony", dims = 1:10) %>%
    FindClusters(resolution = res,algorithm = 1)
}
mpdc <- recluster(subset(seurat_object, celltype == "MP/DC"), res = 0.5)
seurat_object <- drop_cells(seurat_object, WhichCells(mpdc, idents = c("12", "14")))
fib <- recluster(subset(seurat_object, celltype == "Fibroblast"), res = 0.3)
seurat_object <- drop_cells(seurat_object, WhichCells(fib, idents = "6"))
rm(mpdc, fib)
log_step(seurat_object, "Final")

# -----------------------------------------------------------------------------
# 9. Save
# -----------------------------------------------------------------------------
meta <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))
seurat_object$group      <- factor(seurat_object$group, levels = unique(meta$Group))
seurat_object$orig.ident <- factor(seurat_object$orig.ident, levels = meta$orig.ident)

saveRDS(seurat_object, file.path(DATA_DIR, "seurat_object_QC4.rds"))
write.csv(steps, file.path(REPO_DIR, "logs", "cleanup", "cell_filtering_steps.csv"), row.names = FALSE)
write.csv(data.frame(barcode    = colnames(seurat_object),
                     orig.ident = seurat_object$orig.ident,
                     group      = seurat_object$group,
                     celltype   = seurat_object$celltype),
          gzfile(file.path(REPO_DIR, "metadata", "cell_annotations.csv.gz")), row.names = FALSE)

print(steps)
print(table(seurat_object$celltype, seurat_object$group))
sessionInfo()
