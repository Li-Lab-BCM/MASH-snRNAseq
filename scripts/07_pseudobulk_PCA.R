# =============================================================================
# 07_pseudobulk_PCA.R
# Sample-level pseudo-bulk PCA for two comparisons:
#   (1) Chow vs MASH
#   (2) Spp1 f/f (MASH_flox) vs Spp1 f/f;Alb-Cre (MASH_SPP1KO)
#
# Input : <DATA_DIR>/seurat_final.rds
# Output: <DATA_DIR>/results/figures/pseudobulk_PCA_*.pdf
# =============================================================================


library(Seurat)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)


DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")
FIG_DIR  <- file.path(DATA_DIR, "results", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

seurat_object <- readRDS(file.path(DATA_DIR, "seurat_final.rds"))
meta <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))

# Genes excluded: mitochondrial, ribosomal, hemoglobin, immediate-early and sex genes
genes_all <- rownames(seurat_object)
blacklist <- unique(c(
  grep("^mt-", genes_all, value = TRUE),
  grep("^(Rpl|Rps|Mrpl|Mrps)", genes_all, value = TRUE),
  grep("^(Hba|Hbb|Alas2)", genes_all, value = TRUE),
  c("Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3", "Nr4a1", "Nr4a2",
    "Nr4a3", "Btg2", "Dusp1", "Ier2", "Ier3", "Atf3", "Zfp36", "Ccn1", "Nfkbiz"),
  c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty")
))

# Pseudo-bulk counts per sample (log2 CPM)
bulk <- AggregateExpression(seurat_object, group.by = "orig.ident", assays = "RNA",
                            slot = "counts", return.seurat = FALSE)$RNA
colnames(bulk) <- gsub("-", "_", colnames(bulk))
bulk <- bulk[!rownames(bulk) %in% blacklist, ]

# -----------------------------------------------------------------------------
# PCA for a subset of samples; batch correction (limma::removeBatchEffect,
# group kept in the design) is applied only when the samples span >1 batch
# -----------------------------------------------------------------------------
run_pca <- function(samples, labels, group_labels, colors, ellipse = TRUE, n_hvg = 2000) {
  counts <- bulk[, samples]
  counts <- counts[rowSums(counts) > 0, ]
  logcpm <- log2(t(t(counts) / colSums(counts) * 1e6) + 1)

  m <- meta[match(samples, meta$orig.ident), ]
  group <- factor(group_labels[m$Group], levels = unique(group_labels))
  batch <- factor(m$Batch)
  if (nlevels(batch) > 1) {
    logcpm <- removeBatchEffect(logcpm, batch = batch, design = model.matrix(~ group))
  }

  hvg <- names(sort(apply(logcpm, 1, var), decreasing = TRUE))[1:n_hvg]
  pca <- prcomp(scale(t(logcpm[hvg, ])), center = TRUE, scale. = FALSE)
  var_pct <- round(summary(pca)$importance[2, 1:2] * 100, 1)

  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   label = labels[samples], group = group)

  p <- ggplot(df, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey80") +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey80")
  if (ellipse) {
    p <- p + stat_ellipse(aes(colour = group), type = "norm", level = 0.95,
                          linewidth = 0.5, show.legend = FALSE)
  }
  p +
    geom_point(aes(fill = group), shape = 21, size = 2.5, colour = "black", stroke = 0.3) +
    geom_text_repel(aes(label = label), parse = TRUE, size = 4.5,
                    box.padding = 0.6, max.overlaps = Inf) +
    scale_fill_manual(values = colors, labels = parse(text = names(colors)), name = NULL) +
    scale_colour_manual(values = colors) +
    labs(title = "PCA", x = paste0("PC1 (", var_pct[1], "%)"),
         y = paste0("PC2 (", var_pct[2], "%)")) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.border = element_rect(fill = NA, colour = "black"),
          legend.text.align = 0)
}

# (1) Chow vs MASH ------------------------------------------------------------
p1 <- run_pca(
  samples = c("CHOW_1", "CHOW_2", "CHOW_3", "MASH_1", "MASH_2", "MASH_3"),
  labels  = c(CHOW_1 = "Chow1", CHOW_2 = "Chow2", CHOW_3 = "Chow3",
              MASH_1 = "MASH1", MASH_2 = "MASH2", MASH_3 = "MASH3"),
  group_labels = c(CHOW = "Chow", MASH = "MASH"),
  colors  = c(Chow = "#BDBDBD", MASH = "#F4A3A0")
)
ggsave(file.path(FIG_DIR, "pseudobulk_PCA_Chow_vs_MASH.pdf"), p1, width = 4.5, height = 3.5)

# (2) Spp1 f/f vs Spp1 f/f;Alb-Cre ---------------------------------------------
flox <- "italic(Spp1)^{f/f}"
ko   <- "italic(Spp1)^{f/f}*';'*italic(Alb)^{Cre}"
p2 <- run_pca(
  samples = c("MASH_flox_1", "MASH_flox_2", "MASH_SPP1KO_1", "MASH_SPP1KO_2"),
  labels  = c(MASH_flox_1   = paste0(flox, "*1"), MASH_flox_2   = paste0(flox, "*2"),
              MASH_SPP1KO_1 = paste0(ko, "*1"),   MASH_SPP1KO_2 = paste0(ko, "*2")),
  group_labels = setNames(c(flox, ko), c("MASH_flox", "MASH_SPP1KO")),
  colors  = setNames(c("#C2477A", "#4CAF50"), c(flox, ko)),
  ellipse = FALSE
)
ggsave(file.path(FIG_DIR, "pseudobulk_PCA_flox_vs_SPP1KO.pdf"), p2, width = 5, height = 3.5)

sessionInfo()
