# =============================================================================
# Sub-clustering of MP/DC cells (Harmony), sub-type annotation, UMAP and marker plots,
# and export of metadata for RNA velocity (09b).
#
# Input : <DATA_DIR>/seurat_final.rds
# Output: <DATA_DIR>/MPDC_subclusters.rds
#         <DATA_DIR>/results/figures/MPDC_subclusters_umap.pdf
#         <DATA_DIR>/velocity/metadata_MP_velocity.csv
#         metadata/mpdc_subtypes.csv.gz
# =============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)


DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")
FIG_DIR  <- file.path(DATA_DIR, "results", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DATA_DIR, "velocity"), showWarnings = FALSE)

# 1. Sub-cluster MP/DC on the Harmony embedding --------------------------------
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_final.rds"))
mp <- subset(seurat_object, subset = celltype == "MP/DC") %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>%
  FindClusters(resolution = 0.3, algorithm = 1) %>%
  RunUMAP(reduction = "harmony", dims = 1:10, min.dist = 0.5)
rm(seurat_object)

markers <- FindAllMarkers(mp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.4)
write.csv(markers, file.path(FIG_DIR, "MPDC_subcluster_markers.csv"))

# 2. Sub-type annotation ---------------------------------------------------------
#    clusters 3 and 4 share the same resident-macrophage profile and are merged;
#    cluster 5 is split into DCs (Flt3/Xcr1/Clec9a/Zbtb46+, Ccr2-low) and monocytes
cl <- as.character(mp$seurat_clusters)
cl[cl == "4"] <- "3"

c5 <- colnames(mp)[cl == "5"]
expr <- GetAssayData(mp, layer = "data")[c("Flt3", "Xcr1", "Clec9a", "Zbtb46", "Ccr2"), c5]
is_dc <- (expr["Flt3", ] > 0 | expr["Xcr1", ] > 0 | expr["Clec9a", ] > 0 |
          expr["Zbtb46", ] > 0) & expr["Ccr2", ] < 2
cl[cl == "5"] <- ifelse(is_dc, "DC", "Mono")

subtype_map <- c("0" = "MP1", "3" = "MP2", "2" = "MP3", "1" = "MP4", Mono = "Mono", DC = "DC")
mp$celltype <- factor(subtype_map[cl], levels = c("MP1", "MP2", "MP3", "MP4", "Mono", "DC"))
Idents(mp) <- "celltype"
saveRDS(mp, file.path(DATA_DIR, "MPDC_subclusters.rds"))

# 3. UMAP and marker dot plot --------------------------------------------------
mp_colors <- c(MP1 = "#E7298A", MP2 = "#66A61E", MP3 = "#E6AB02",
               MP4 = "#1B9E77", Mono = "#7570B3", DC = "#D95F02")
p <- DimPlot(mp, reduction = "umap", cols = mp_colors, label = TRUE,
             repel = TRUE, pt.size = 0.5) +
  ggtitle("MP/DC subclusters") + theme(plot.title = element_text(hjust = 0.5))
ggsave(file.path(FIG_DIR, "MPDC_subclusters_umap.pdf"), p, width = 8, height = 6)

marker_list <- c(
  "Ccr2", "H2-Ab1", "H2-Aa", "Tgfbr1", "Ciita",
  "Vsig4", "Lyz2", "C1qa", "C1qc", "Cd163",
  "Mrc1", "Adgre1", "Cd68", "Mertk",
  "Timd4", "Lyve1", "Folr2", "Cbr2", "Pf4",
  "Ly6c2", "Plac8", "Cx3cr1", "Fcn1", "Cd244a",
  "Flt3", "Xcr1", "Clec9a", "Irf8", "Zbtb46"
)
p_dot <- DotPlot(mp, features = intersect(marker_list, rownames(mp)), dot.scale = 6) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, name = "Scaled Average\nExpression") +
  scale_y_discrete(limits = rev(levels(mp$celltype))) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 10))
ggsave(file.path(FIG_DIR, "MPDC_marker_dotplot.pdf"), p_dot, width = 7.5, height = 3)


# 3b. Figure 3F, Figure 3I, and anti-OPN MP1 inflammatory-gene DotPlots --------
# Dot color = scaled average expression (Seurat DotPlot, scale = TRUE)
# Dot size  = percentage of nuclei expressing the indicated gene

make_pair_dotplot <- function(object, genes, groups, plot_title = NULL) {
  x <- subset(object, subset = group %in% groups)
  x$group <- factor(x$group, levels = groups)

  DotPlot(x, features = genes, group.by = "group", scale = TRUE) +
    scale_colour_gradient2(
      low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
      midpoint = 0, name = "Scaled Average\nExpression"
    ) +
    scale_size_continuous(name = "Pct. expressed", range = c(2, 9)) +
    geom_point(
      aes(size = pct.exp),
      shape = 21, colour = "black", fill = NA, stroke = 1
    ) +
    coord_flip() +
    ggtitle(plot_title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
      panel.grid = element_blank()
    )
}

# Figure 3F: whole MP/DC cluster
genes_fig3f <- c("Ccr2", "Il1b", "Casp1", "Tgfbr1", "Pycard", "Aim2", "Trem2")

p3f_disease <- make_pair_dotplot(
  mp, genes_fig3f, c("CHOW", "MASH"), "CHOW vs MASH"
)
p3f_ko <- make_pair_dotplot(
  mp, genes_fig3f, c("MASH_flox", "MASH_SPP1KO"),
  "MASH_flox vs MASH_SPP1KO"
)
p3f <- p3f_disease + p3f_ko + plot_layout(ncol = 2)
ggsave(
  file.path(FIG_DIR, "Figure3F_MPDC_inflammatory_dotplot.pdf"),
  p3f, width = 8.5, height = 5
)

# Figure 3I: MP1 subcluster
mp1 <- subset(mp, subset = celltype == "MP1")
genes_fig3i <- c("Trem2", "Tgfbr1", "Pycard", "Aim2", "Tlr1", "Traf6")

p3i_disease <- make_pair_dotplot(
  mp1, genes_fig3i, c("CHOW", "MASH"), "CHOW vs MASH"
)
p3i_ko <- make_pair_dotplot(
  mp1, genes_fig3i, c("MASH_flox", "MASH_SPP1KO"),
  "MASH_flox vs MASH_SPP1KO"
)
p3i <- p3i_disease + p3i_ko + plot_layout(ncol = 2)
ggsave(
  file.path(FIG_DIR, "Figure3I_MP1_inflammatory_dotplot.pdf"),
  p3i, width = 8.5, height = 5
)

# Figure 7L / anti-OPN MP1 comparison
p_mp1_anti <- make_pair_dotplot(
  mp1, genes_fig3i, c("MASH", "MASHantiSPP1"),
  "MASHantiSPP1 vs MASH"
)
ggsave(
  file.path(FIG_DIR, "MP1_MASHantiSPP1_vs_MASH_dotplot.pdf"),
  p_mp1_anti, width = 4, height = 5
)

# 3c. Supplemental Figure S7H: MP/DC subtype proportions per sequenced sample --
sample_order <- c(
  "CHOW_1", "CHOW_2", "CHOW_3",
  "MASH_1", "MASH_2", "MASH_3",
  "MASH_flox_1", "MASH_flox_2",
  "MASH_SPP1KO_1", "MASH_SPP1KO_2"
)

prop_df <- mp@meta.data %>%
  filter(orig.ident %in% sample_order) %>%
  mutate(orig.ident = factor(as.character(orig.ident), levels = sample_order)) %>%
  count(orig.ident, celltype, name = "n") %>%
  group_by(orig.ident) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

p_s7h <- ggplot(prop_df, aes(orig.ident, proportion, fill = celltype)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  scale_fill_manual(values = mp_colors, name = "MP/DC subtype") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Proportion of MP/DC nuclei") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(
  file.path(FIG_DIR, "Supplement_Figure_S7H_MPDC_proportion_by_sample.pdf"),
  p_s7h, width = 8.5, height = 5.2
)

# 4. Exports -------------------------------------------------------------------
write.csv(data.frame(barcode = colnames(mp), orig.ident = mp$orig.ident,
                     group = mp$group, subtype = mp$celltype),
          gzfile(file.path(REPO_DIR, "metadata", "mpdc_subtypes.csv.gz")), row.names = FALSE)

vel <- subset(mp, subset = celltype != "DC")          # DCs excluded from velocity
harmony <- Embeddings(vel, "harmony")
colnames(harmony) <- paste0("harmony_", seq_len(ncol(harmony)))
umap <- Embeddings(vel, "umap")
colnames(umap) <- c("UMAP_1", "UMAP_2")
write.csv(data.frame(cell_id = colnames(vel), orig.ident = vel$orig.ident,
                     group = vel$group, celltype = vel$celltype,
                     umap, harmony, check.names = FALSE),
          file.path(DATA_DIR, "velocity", "metadata_MP_velocity.csv"), row.names = FALSE)

print(table(mp$celltype, mp$group))
sessionInfo()
