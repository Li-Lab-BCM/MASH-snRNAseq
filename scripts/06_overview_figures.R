# =============================================================================
# 06_overview_figures.R
# Final Harmony embedding of the annotated dataset, cell-type UMAP,
# cell-type proportions and marker dot plot.
#
# Input : <DATA_DIR>/seurat_object_QC4.rds
# Output: <DATA_DIR>/seurat_final.rds
#         <DATA_DIR>/results/figures/
# =============================================================================


library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

DATA_DIR <- "/path/to/MASHproject"
FIG_DIR  <- file.path(DATA_DIR, "results", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1234)

cell_type_colors <- c(
  "Endocardial cell"     = "#AE3EC9FF", "Epicardial cell" = "#4DBBD5FF",
  "Fibroblast"           = "#EFC000FF", "MP/DC"           = "#DC0000FF",
  "Endothelial cell"     = "#3C5488FF", "Atrial Cardiomyocyte" = "#8491B4FF",
  "Mural cell"           = "#91D1C2FF", "T/NK cell"       = "#F39B7FFF",
  "LEC"                  = "#7E6148FF", "Adipocyte"       = "#808180FF",
  "Schwann cell"         = "pink",      "B cell"          = "#00A087FF"
)

# 1. Harmony integration of the final dataset ---------------------------------
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_object_QC4.rds"))
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object) %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunHarmony("orig.ident")
seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, min.dist = 0.5)
saveRDS(seurat_object, file.path(DATA_DIR, "seurat_final.rds"))

# 2. Cell-type UMAP ------------------------------------------------------------
p <- DimPlot(seurat_object, reduction = "umap", group.by = "celltype",
             pt.size = 0.3, raster = FALSE, cols = cell_type_colors) + theme_void()
ggsave(file.path(FIG_DIR, "UMAP_celltype.pdf"), p, width = 8, height = 6)

# 3. Cell-type proportions per sample and per group -------------------------
plot_proportion <- function(obj, x_col, xlab, x_levels = NULL, x_labels = NULL) {
  df <- obj@meta.data %>%
    transmute(
      celltype = as.character(celltype),
      x = as.character(.data[[x_col]])
    ) %>%
    count(celltype, x, name = "n")

  if (!is.null(x_levels)) {
    df$x <- factor(df$x, levels = x_levels)
  }

  p <- ggplot(df, aes(x = x, y = n, fill = celltype)) +
    geom_col(position = "fill", width = 0.7) +
    scale_fill_manual(values = cell_type_colors) +
    labs(x = xlab, y = "Proportion", fill = "Cell type") +
    theme_minimal(base_size = 15) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(x_labels)) {
    p <- p + scale_x_discrete(labels = x_labels)
  }

  p
}

ggsave(
  file.path(FIG_DIR, "proportion_by_sample.pdf"),
  plot_proportion(seurat_object, "orig.ident", "Sample"),
  width = 8,
  height = 6
)

ggsave(
  file.path(FIG_DIR, "proportion_by_group.pdf"),
  plot_proportion(seurat_object, "group", "Group"),
  width = 5,
  height = 6
)

# Supplemental Figure S13A: Chow, MASH and anti-OPN ---------------------------
s13a <- subset(
  seurat_object,
  subset = group %in% c("CHOW", "MASH", "MASHantiSPP1")
)

s13a$group <- factor(
  s13a$group,
  levels = c("CHOW", "MASH", "MASHantiSPP1")
)

p_s13a_umap <- DimPlot(
  s13a,
  reduction = "umap",
  group.by = "celltype",
  split.by = "group",
  pt.size = 0.25,
  raster = FALSE,
  cols = cell_type_colors
) +
  theme_void()

ggsave(
  file.path(FIG_DIR, "Supplement_Figure_S13A_UMAP.pdf"),
  p_s13a_umap,
  width = 15,
  height = 5
)

p_s13a_prop <- plot_proportion(
  s13a,
  "group",
  "Group",
  x_levels = c("CHOW", "MASH", "MASHantiSPP1"),
  x_labels = c(
    CHOW = "Chow",
    MASH = "MASH",
    MASHantiSPP1 = "MASH + anti-OPN"
  )
)

ggsave(
  file.path(FIG_DIR, "Supplement_Figure_S13A_celltype_proportion.pdf"),
  p_s13a_prop,
  width = 5.5,
  height = 6
)

# 4. Marker dot plot -----------------------------------------------------------
markers <- c(
  "Car3", "Adipoq", "Cidec",        # Adipocyte
  "Tnnt2", "Nppa", "Myl4",          # Atrial cardiomyocyte
  "Cd79b", "Cd79a", "Cd19",         # B cell
  "H19", "Erg", "Npr3",             # Endocardial cell
  "Pecam1", "Tek", "Cdh5",          # Endothelial cell
  "Tbx18", "Wt1", "Upk1b",          # Epicardial cell
  "Col1a1", "Gsn", "Pdgfra",        # Fibroblast
  "Prox1", "Mmrn1", "Reln",         # LEC
  "C1qb", "C1qc", "Itgam",          # MP/DC
  "Pdgfrb", "Rgs5", "Cspg4",        # Mural cell
  "Plp1", "Sox10", "Lgi4",          # Schwann cell
  "Cd3e", "Cd3g", "Cd4"             # T/NK cell
)
celltype_order <- c("Adipocyte", "Atrial Cardiomyocyte", "B cell", "Endocardial cell",
                    "Endothelial cell", "Epicardial cell", "Fibroblast", "LEC",
                    "MP/DC", "Mural cell", "Schwann cell", "T/NK cell")
seurat_object$celltype <- factor(seurat_object$celltype, levels = celltype_order)

p <- DotPlot(seurat_object, features = intersect(markers, rownames(seurat_object)),
             group.by = "celltype", dot.scale = 6) +
  coord_flip() + RotatedAxis() +
  scale_color_gradientn(colors = c("lightblue", "white", "lightyellow", "red"))
ggsave(file.path(FIG_DIR, "marker_dotplot.pdf"), p, width = 6, height = 9)

sessionInfo()
