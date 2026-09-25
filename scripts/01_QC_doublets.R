# =============================================================================
# 01_QC_doublets.R
# Per-sample QC filtering and doublet removal (scDblFinder) on CellBender output.
#
# Input : <DATA_DIR>/<Sample>/cellbender_cleaned_filtered.h5ad
# Output: <DATA_DIR>/<Sample>-cellbender-SCDBL_final.rds
#         <DATA_DIR>/QC_plots/<Sample>_*.pdf
#         logs/qc/QC_summary_all_samples.csv
# =============================================================================


library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
library(zellkonverter)
library(patchwork)
library(dplyr)
library(ggplot2)

DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")
PLOT_DIR <- file.path(DATA_DIR, "QC_plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(REPO_DIR, "logs", "qc"), showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))

qc_params <- list(
  nFeature_min = 500,
  nFeature_max = 6000,
  nCount_min   = 1000,
  nCount_max   = 40000,
  mt_max       = 10,   # % mitochondrial
  hb_max       = 1,    # % hemoglobin (column name: percent.rb)
  dbl_warn     = 15    # warn if doublet rate (%) exceeds this
)

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "percent.rb")
qc_summary  <- list()

for (sample in meta$Sample) {
  message("Processing ", sample)
  set.seed(1234)

  # 1. Load CellBender output --------------------------------------------------
  sce_raw <- readH5AD(file.path(DATA_DIR, sample, "cellbender_cleaned_filtered.h5ad"))
  A <- CreateSeuratObject(counts = assays(sce_raw)[["X"]], project = sample,
                          min.cells = 3, min.features = 200)
  n_raw <- ncol(A)

  # 2. QC metrics --------------------------------------------------------------
  A[["percent.mt"]]   <- PercentageFeatureSet(A, pattern = "^mt-")
  A[["percent.ribo"]] <- PercentageFeatureSet(A, pattern = "^Rps|^Rpl")
  A[["percent.rb"]]   <- PercentageFeatureSet(A, pattern = "^Hba|^Hbb")  # hemoglobin

  pdf(file.path(PLOT_DIR, paste0(sample, "_QC_violin_prefilter.pdf")), width = 14, height = 5)
  print(VlnPlot(A, features = qc_features, ncol = 5, pt.size = 0) +
          plot_annotation(title = paste0(sample, " - Pre-filter QC")))
  dev.off()

  pdf(file.path(PLOT_DIR, paste0(sample, "_QC_scatter.pdf")), width = 7, height = 6)
  print(FeatureScatter(A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
          geom_hline(yintercept = c(qc_params$nFeature_min, qc_params$nFeature_max),
                     linetype = "dashed", color = "red") +
          geom_vline(xintercept = c(qc_params$nCount_min, qc_params$nCount_max),
                     linetype = "dashed", color = "red") +
          ggtitle(paste0(sample, " - nCount vs nFeature")))
  dev.off()

  # 3. QC filtering (before doublet detection) ---------------------------------
  A <- subset(A, subset =
                nFeature_RNA > qc_params$nFeature_min &
                nFeature_RNA < qc_params$nFeature_max &
                nCount_RNA   > qc_params$nCount_min   &
                nCount_RNA   < qc_params$nCount_max   &
                percent.mt   < qc_params$mt_max       &
                percent.rb   < qc_params$hb_max)
  n_post_qc <- ncol(A)

  # 4. Preprocessing for doublet detection / visualisation ---------------------
  A <- NormalizeData(A) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30, verbose = FALSE)

  # 5. Doublet detection -------------------------------------------------------
  sce <- scDblFinder(as.SingleCellExperiment(A))
  A$doublet_status <- sce$scDblFinder.class
  A$doublet_score  <- sce$scDblFinder.score
  dbl_rate <- round(mean(A$doublet_status == "doublet") * 100, 1)
  if (dbl_rate > qc_params$dbl_warn) warning(sample, ": high doublet rate (", dbl_rate, "%)")

  A <- RunUMAP(A, dims = 1:20)
  pdf(file.path(PLOT_DIR, paste0(sample, "_UMAP_doublet.pdf")), width = 16, height = 5)
  print(
    DimPlot(A, group.by = "doublet_status",
            cols = c(singlet = "steelblue", doublet = "red")) +
      ggtitle(paste0(sample, ": Doublets (", dbl_rate, "%)")) |
    FeaturePlot(A, features = "doublet_score") |
    FeaturePlot(A, features = "percent.mt")
  )
  dev.off()

  # 6. Remove doublets ---------------------------------------------------------
  A <- subset(A, subset = doublet_status == "singlet")
  n_final <- ncol(A)

  pdf(file.path(PLOT_DIR, paste0(sample, "_QC_violin_postfilter.pdf")), width = 14, height = 5)
  print(VlnPlot(A, features = qc_features, ncol = 5, pt.size = 0) +
          plot_annotation(title = paste0(sample, " - Post-filter QC")))
  dev.off()

  saveRDS(A, file.path(DATA_DIR, paste0(sample, "-cellbender-SCDBL_final.rds")))

  qc_summary[[sample]] <- data.frame(
    sample          = sample,
    n_raw           = n_raw,
    n_post_qc       = n_post_qc,
    n_doublet       = n_post_qc - n_final,
    dbl_rate_pct    = dbl_rate,
    n_final         = n_final,
    median_nFeature = median(A$nFeature_RNA),
    median_nCount   = median(A$nCount_RNA),
    median_mt       = round(median(A$percent.mt), 2)
  )
  message(sample, ": ", n_raw, " -> QC ", n_post_qc, " -> final ", n_final,
          " (doublets ", dbl_rate, "%)")
}

qc_df <- do.call(rbind, qc_summary)
write.csv(qc_df, file.path(REPO_DIR, "logs", "qc", "QC_summary_all_samples.csv"), row.names = FALSE)
print(qc_df)

sessionInfo()
