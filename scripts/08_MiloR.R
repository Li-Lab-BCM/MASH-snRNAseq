# =============================================================================
# Differential abundance of cell neighbourhoods (miloR) across all cell types.
# Design: ~ batch + group; contrasts MASH vs CHOW, MASHantiSPP1 vs MASH,
# MASH_SPP1KO vs MASH_flox.
#
# Input : <DATA_DIR>/seurat_final.rds
# Output: <DATA_DIR>/results/MiloR/
# =============================================================================


library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)

DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")
OUT_DIR  <- file.path(DATA_DIR, "results", "MiloR")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1234)   # neighbourhood sampling in makeNhoods is random

seurat_object <- readRDS(file.path(DATA_DIR, "seurat_final.rds"))
meta <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))

# 1. Milo object, kNN graph and neighbourhoods on the Harmony embedding --------
milo <- Milo(as.SingleCellExperiment(JoinLayers(seurat_object)))
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = "HARMONY")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, reduced_dims = "HARMONY")
ggsave(file.path(OUT_DIR, "nhood_size_hist.pdf"), plotNhoodSizeHist(milo), width = 6, height = 4)

milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = "orig.ident")
milo <- calcNhoodDistance(milo, d = 30, reduced.dim = "HARMONY")

# 2. Design -------------------------------------------------------------------
sample_info <- data.frame(
  sample = meta$orig.ident,
  group  = factor(meta$Group, levels = unique(meta$Group)),
  batch  = factor(meta$Batch),
  row.names = meta$orig.ident
)[colnames(nhoodCounts(milo)), ]
design <- model.matrix(~ batch + group, data = sample_info)

contrasts <- c(
  MASH_vs_CHOW             = "groupMASH",
  MASHantiSPP1_vs_MASH     = "groupMASHantiSPP1 - groupMASH",
  MASH_SPP1KO_vs_MASH_flox = "groupMASH_SPP1KO - groupMASH_flox"
)

# 3. Differential abundance testing -------------------------------------------
milo <- buildNhoodGraph(milo)
da_list <- list()
for (nm in names(contrasts)) {
  da <- testNhoods(milo, design = design, design.df = sample_info,
                   reduced.dim = "HARMONY", model.contrasts = contrasts[[nm]],
                   fdr.weighting = "graph-overlap")
  da <- annotateNhoods(milo, da, coldata_col = "celltype")
  da_list[[nm]] <- da
  write.csv(da, file.path(OUT_DIR, paste0("DA_", nm, ".csv")), row.names = FALSE)

  ggsave(file.path(OUT_DIR, paste0("DA_beeswarm_", nm, ".pdf")),
         plotDAbeeswarm(da, group.by = "celltype", alpha = 1) + ggtitle(nm) +
           scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                 midpoint = 0, limits = c(-4, 4), oob = scales::squish,
                                 name = "logFC") +
           theme(axis.text.y = element_text(size = 8)),
         width = 8, height = 6)
  ggsave(file.path(OUT_DIR, paste0("DA_umap_", nm, ".pdf")),
         plotNhoodGraphDA(milo, da, layout = "UMAP", alpha = 1) +
           scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                midpoint = 0, limits = c(-4, 4), oob = scales::squish,
                                name = "logFC") +
           ggtitle(nm),
         width = 8, height = 6)
}

# 4. Summary of significant neighbourhoods (SpatialFDR < 0.1) by cell type --
da_summary <- bind_rows(lapply(names(da_list), function(nm) {
  da_list[[nm]] %>%
    filter(!is.na(SpatialFDR), SpatialFDR < 0.1) %>%
    mutate(direction = ifelse(logFC > 0, "enriched", "depleted")) %>%
    group_by(celltype, direction) %>%
    summarise(n_nhoods = n(), mean_logFC = mean(logFC), .groups = "drop") %>%
    mutate(comparison = nm)
}))
write.csv(da_summary, file.path(OUT_DIR, "DA_summary_by_celltype.csv"), row.names = FALSE)
print(sapply(da_list, function(x) sum(x$SpatialFDR < 0.1, na.rm = TRUE)))
saveRDS(da_list, file.path(OUT_DIR, "da_results_list.rds"))
saveRDS(milo, file.path(OUT_DIR, "milo_object.rds"))

sessionInfo()
