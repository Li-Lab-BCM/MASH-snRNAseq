# =============================================================================
# Fibroblast analyses: sub-clusters and markers, activated-fibroblast (AFB)
# score, sample-level module scores, and pseudo-bulk DESeq2 with quadrant plot.
#
# Input : <DATA_DIR>/seurat_final.rds
# Output: <DATA_DIR>/FB_subclusters.rds
#         <DATA_DIR>/results/fibroblast/
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)
library(msigdbr)
library(Nebulosa)
library(Matrix)
library(edgeR)
library(DESeq2)
library(ggrepel)

DATA_DIR <- "/path/to/MASHproject"
REPO_DIR <- file.path(DATA_DIR, "MASH-snRNAseq")
OUT_DIR  <- file.path(DATA_DIR, "results", "fibroblast")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1234)

group_order <- c("CHOW", "MASH", "MASHantiSPP1", "MASH_flox", "MASH_SPP1KO")
fb_colors <- c(FB1 = "#FF7F00", FB2 = "#A65628", FB3 = "#4DAF4A",
               FB4 = "#984EA3", FB5 = "#377EB8")
group_colors <- c(CHOW = "#BFBFBF", MASH = "#D73027", MASHantiSPP1 = "#CCA961",
                  MASH_flox = "#C95D7B", MASH_SPP1KO = "#1A9850")

# Axis labels used in the figures
group_labels <- c(CHOW = "Chow", MASH = "MASH", MASHantiSPP1 = "anti-OPN",
                  MASH_flox = "italic(Spp1)^{f/f}",
                  MASH_SPP1KO = "italic(Spp1)^{f/f}*';'*italic(Alb)^{Cre}")
label_groups <- function(x) parse(text = group_labels[x])
fig_groups <- c("CHOW", "MASH", "MASH_flox", "MASH_SPP1KO")   # groups shown in the figure

# =============================================================================
# 1. Sub-clustering (Harmony) and naming
# =============================================================================
seurat_object <- readRDS(file.path(DATA_DIR, "seurat_final.rds"))
FB <- subset(seurat_object, subset = celltype == "Fibroblast") %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.3,algorithm = 1) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, min.dist = 0.3)
rm(seurat_object)

write.csv(FindAllMarkers(FB, logfc.threshold = 0.4, min.pct = 0.1, only.pos = TRUE),
          file.path(OUT_DIR, "FB_subcluster_markers.csv"))

fb_map <- c("1" = "FB1", "0" = "FB2", "2" = "FB3", "3" = "FB4", "4" = "FB5")
FB$fb_cluster <- factor(fb_map[as.character(FB$seurat_clusters)], levels = names(fb_colors))
FB$group <- factor(FB$group, levels = group_order)
Idents(FB) <- "fb_cluster"

p <- DimPlot(FB, reduction = "umap", label = TRUE, pt.size = 0.5, cols = fb_colors) +
  ggtitle("Fibroblast subclusters")
ggsave(file.path(OUT_DIR, "FB_subclusters_umap.pdf"), p, width = 6, height = 6)

# Marker dot plot: 3 representative markers per sub-cluster
# (Chow, MASH, Spp1 f/f and Spp1 f/f;Alb-Cre nuclei)
marker_list <- list(FB1 = c("Id1", "Dkk3", "Piezo2"),
                    FB2 = c("Nrxn1", "Kcnc2", "Notch3"),
                    FB3 = c("Nkain2", "Lvrn", "Ccbe1"),
                    FB4 = c("Limch1", "Pi16", "Hmcn2"),
                    FB5 = c("Ch25h", "Nr4a3", "Thbs1"))
markers_use <- intersect(unlist(marker_list, use.names = FALSE), rownames(FB))
fb_marker <- subset(FB, subset = group %in% fig_groups)
dp <- DotPlot(fb_marker, features = markers_use, group.by = "fb_cluster", dot.scale = 7)$data
dp$features.plot <- factor(dp$features.plot, levels = markers_use)
dp$id <- factor(dp$id, levels = rev(names(fb_colors)))
p <- ggplot(dp, aes(features.plot, id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5), linewidth = 0.35, color = "grey80") +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                        name = "Scaled average\nexpression") +
  scale_size_continuous(range = c(1.5, 8), name = "Percent expressed") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", color = "black"),
        axis.text.y = element_text(color = "black"), axis.ticks = element_blank())
ggsave(file.path(OUT_DIR, "FB_marker_dotplot.pdf"), p, width = 8.5, height = 4)
rm(fb_marker)

# Sub-cluster proportions per group
plot_fb_proportion <- function(groups, labels = TRUE) {
  df <- FB@meta.data %>%
    filter(group %in% groups) %>%
    mutate(group = factor(group, levels = groups)) %>%
    count(group, fb_cluster) %>%
    group_by(group) %>%
    mutate(prop = n / sum(n),
           label = ifelse(prop >= 0.05, paste0(round(prop * 100, 1), "%"), "")) %>%
    ungroup()
  p <- ggplot(df, aes(group, prop, fill = fb_cluster)) +
    geom_col(width = 0.8, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = fb_colors, name = "Fibroblast\nsubcluster") +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(title = "Proportion of Fibroblast Subclusters", x = NULL, y = "Proportion") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))
  if (labels) p <- p + geom_text(aes(label = label), position = position_stack(vjust = 0.5),
                                 size = 3.5, color = "white", fontface = "bold")
  p
}
ggsave(file.path(OUT_DIR, "FB_proportion_5groups.pdf"),
       plot_fb_proportion(group_order), width = 4, height = 5.5)
ggsave(file.path(OUT_DIR, "FB_proportion_noAntiSPP1.pdf"),
       plot_fb_proportion(c("CHOW", "MASH", "MASH_flox", "MASH_SPP1KO"), labels = FALSE),
       width = 4, height = 5.2)

# =============================================================================
# 2. Activated fibroblast (AFB) score
# =============================================================================
afb_genes <- c("Postn", "Thbs4", "Cthrc1", "Comp", "Ltbp2", "Col1a1",
               "Timp1", "Adam12", "Acta2", "Tagln", "Meox1")
FB <- AddModuleScore(FB, features = list(intersect(afb_genes, rownames(FB))), name = "AFB")

afb <- data.frame(fb_cluster = FB$fb_cluster, group = FB$group, score = FB$AFB1) %>%
  filter(group %in% fig_groups) %>%
  mutate(group = factor(group, levels = fig_groups),
         score_z = as.numeric(scale(score))) %>%
  group_by(fb_cluster, group) %>%
  summarise(n = n(), mean_score = mean(score), mean_score_z = mean(score_z),
            pct_pos = 100 * mean(score > 0), .groups = "drop")
cl_order <- afb %>% group_by(fb_cluster) %>% summarise(k = mean(mean_score)) %>%
  arrange(desc(k)) %>% pull(fb_cluster)
afb$fb_cluster <- factor(afb$fb_cluster, levels = cl_order)
write.csv(afb, file.path(OUT_DIR, "AFB_score_summary.csv"), row.names = FALSE)

p <- ggplot(afb, aes(group, fb_cluster)) +
  geom_point(aes(size = pct_pos, fill = mean_score_z), shape = 21, colour = "black") +
  scale_size(range = c(2, 10), name = "Pct.\npositive (%)") +
  scale_fill_gradientn(colours = c("#313695", "#74add1", "#F7F7F7", "#f46d43", "#a50026"),
                       limits = c(-1.5, 1.5), oob = squish, name = "Mean\nAFB (z)") +
  scale_x_discrete(labels = label_groups) +
  labs(x = NULL, y = NULL, title = "Fibroblast activation score") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "AFB_score_bubble.pdf"), p, width = 5.5, height = 4.5)

p <- plot_density(FB, features = "AFB1", reduction = "umap", pal = "plasma")
ggsave(file.path(OUT_DIR, "AFB_density.pdf"), p, width = 6, height = 6)

# =============================================================================
# 3. Module scores: per-sample means and between-group comparisons
# =============================================================================
msig <- bind_rows(
  msigdbr(species = "Mus musculus", category = "H"),
  msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME")
)
geneset <- function(name) unique(msig$gene_symbol[msig$gs_name == name])

signatures <- list(
  TGFb_Hallmark              = geneset("HALLMARK_TGF_BETA_SIGNALING"),
  EMT_Hallmark               = geneset("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
  FB_inflammatory_Hallmark   = geneset("HALLMARK_INFLAMMATORY_RESPONSE"),
  FB_TNFA_NFKB_Hallmark      = geneset("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
  FB_IL6_JAK_STAT3_Hallmark  = geneset("HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  IL1_response_MSigDB        = geneset("REACTOME_INTERLEUKIN_1_SIGNALING"),
  Cellular_senescence_MSigDB = geneset("REACTOME_CELLULAR_SENESCENCE"),
  SASP_MSigDB                = geneset("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP")
)
category <- c(TGFb_Hallmark = "Fibrosis", EMT_Hallmark = "Fibrosis",
              FB_inflammatory_Hallmark = "Inflammation", FB_TNFA_NFKB_Hallmark = "Inflammation",
              FB_IL6_JAK_STAT3_Hallmark = "Inflammation", IL1_response_MSigDB = "IL1_pathway",
              Cellular_senescence_MSigDB = "SASP_senescence", SASP_MSigDB = "SASP_senescence")

write.csv(data.frame(signature = names(signatures), category = category[names(signatures)],
                     n_genes = lengths(signatures),
                     genes = sapply(signatures, paste, collapse = ";")),
          file.path(OUT_DIR, "module_score_signatures.csv"), row.names = FALSE)

for (sn in names(signatures)) {
  FB <- AddModuleScore(FB, features = list(intersect(signatures[[sn]], rownames(FB))),
                       name = paste0(sn, "_MS"), ctrl = 100)
}

# Per-sample mean score
sample_scores <- FB@meta.data %>%
  select(group, sample_id = orig.ident, all_of(paste0(names(signatures), "_MS1"))) %>%
  pivot_longer(-c(group, sample_id), names_to = "signature", values_to = "score") %>%
  mutate(signature = sub("_MS1$", "", signature)) %>%
  group_by(group, sample_id, signature) %>%
  summarise(mean_score = mean(score), n_cells = n(), .groups = "drop")
write.csv(sample_scores, file.path(OUT_DIR, "module_score_sample_level.csv"), row.names = FALSE)

# Wilcoxon test on per-sample means, BH-adjusted across the three comparisons
comparisons <- list(Disease = c("CHOW", "MASH"),
                    Ab_rescue = c("MASH", "MASHantiSPP1"),
                    KO_rescue = c("MASH_flox", "MASH_SPP1KO"))
p_label <- function(p) cut(p, c(-Inf, 1e-4, 1e-3, 0.01, 0.05, Inf),
                           labels = c("****", "***", "**", "*", "ns"))
stats <- sample_scores %>%
  group_by(signature) %>%
  group_modify(~ bind_rows(lapply(names(comparisons), function(cp) {
    x <- .x$mean_score[.x$group == comparisons[[cp]][1]]
    y <- .x$mean_score[.x$group == comparisons[[cp]][2]]
    data.frame(comp = cp, delta = mean(y) - mean(x),
               p_value = wilcox.test(x, y, exact = FALSE)$p.value)
  }))) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         p_label = as.character(p_label(p_adj))) %>%
  ungroup()
write.csv(stats, file.path(OUT_DIR, "module_score_stats.csv"), row.names = FALSE)

# Heatmap of group means (z-scored across the groups shown); statistics in module_score_stats.csv
sig_labels <- c(TGFb_Hallmark = "TGFbeta", EMT_Hallmark = "EMT",
                IL1_response_MSigDB = "IL1_response", FB_inflammatory_Hallmark = "FB_inflammation",
                FB_TNFA_NFKB_Hallmark = "FB_TNFalpha_NFKB",
                FB_IL6_JAK_STAT3_Hallmark = "FB_IL6_JAK_STAT3",
                Cellular_senescence_MSigDB = "Cellular_senescence", SASP_MSigDB = "SASP")

plot_score_heatmap <- function(groups, zscore = TRUE, title) {
  heat <- sample_scores %>%
    filter(group %in% groups) %>%
    mutate(group = factor(group, levels = groups)) %>%
    group_by(signature, group) %>%
    summarise(grp_mean = mean(mean_score), .groups = "drop") %>%
    group_by(signature) %>%
    mutate(fill = if (zscore) as.numeric(scale(grp_mean)) else grp_mean) %>%
    ungroup() %>%
    mutate(signature = factor(sig_labels[signature], levels = rev(sig_labels)))
  ggplot(heat, aes(group, signature, fill = fill)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = if (zscore) "z-score\n(group mean)" else "Group mean") +
    scale_x_discrete(labels = label_groups) +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(color = "black"), panel.grid = element_blank())
}

ggsave(file.path(OUT_DIR, "module_score_heatmap.pdf"),
       plot_score_heatmap(fig_groups, title = "Altered functions in FB"),
       width = 6, height = 5.5)
ggsave(file.path(OUT_DIR, "module_score_heatmap_antiSPP1.pdf"),
       plot_score_heatmap(c("MASH", "MASHantiSPP1"), zscore = FALSE,
                          title = "FB module score: anti-OPN vs MASH"),
       width = 5, height = 5.5)

# Sample-level box plots: AFB and IL-1 response, MASH vs MASHantiSPP1
box_df <- FB@meta.data %>%
  select(group, sample_id = orig.ident, AFB1, IL1_response_MSigDB_MS1) %>%
  filter(group %in% c("MASH", "MASHantiSPP1")) %>%
  pivot_longer(c(AFB1, IL1_response_MSigDB_MS1), names_to = "signature", values_to = "score") %>%
  group_by(group, sample_id, signature) %>%
  summarise(mean_score = mean(score), .groups = "drop") %>%
  mutate(group = factor(group, levels = c("MASH", "MASHantiSPP1")),
         signature = factor(signature, levels = c("AFB1", "IL1_response_MSigDB_MS1"),
                            labels = c("Fibroblast activation", "IL-1 response")))
p <- ggplot(box_df, aes(group, mean_score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 2.5) +
  geom_text(aes(label = sample_id), hjust = -0.1, size = 3) +
  scale_fill_manual(values = group_colors) +
  facet_wrap(~ signature, scales = "free_y", nrow = 1) +
  labs(title = "FB module scores: MASHantiSPP1 vs MASH", x = NULL,
       y = "Sample-level mean score") +
  theme_classic(base_size = 12) +
  theme(strip.text = element_text(face = "bold"), legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(OUT_DIR, "module_score_box_antiSPP1.pdf"), p, width = 7, height = 4.5)


# =============================================================================
# Figure 7N — Fibroblast activation score
# CHOW vs MASH vs anti-OPN
# One point represents one independently sequenced sample/library.
# =============================================================================

fig7n_df <- FB@meta.data %>%
  transmute(
    group = as.character(group),
    sample_id = as.character(orig.ident),
    score = AFB1
  ) %>%
  filter(group %in% c("CHOW", "MASH", "MASHantiSPP1")) %>%
  group_by(group, sample_id) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    group = factor(
      group,
      levels = c("CHOW", "MASH", "MASHantiSPP1")
    )
  )

fig7n_comparisons <- list(
  c("CHOW", "MASH"),
  c("MASH", "MASHantiSPP1")
)

fig7n_stats <- bind_rows(lapply(fig7n_comparisons, function(cmp) {
  x <- fig7n_df %>% filter(group == cmp[1]) %>% pull(mean_score)
  y <- fig7n_df %>% filter(group == cmp[2]) %>% pull(mean_score)

  p_value <- if (length(x) >= 2 && length(y) >= 2) {
    wilcox.test(x, y, exact = FALSE)$p.value
  } else {
    NA_real_
  }

  data.frame(
    group1 = cmp[1],
    group2 = cmp[2],
    p_value = p_value,
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    p_label = case_when(
      is.na(p_adj) ~ "NA",
      p_adj < 0.0001 ~ "****",
      p_adj < 0.001 ~ "***",
      p_adj < 0.01 ~ "**",
      p_adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

fig7n_y_range <- range(fig7n_df$mean_score, na.rm = TRUE)
fig7n_y_span <- diff(fig7n_y_range)
if (fig7n_y_span == 0) fig7n_y_span <- 1

fig7n_stats <- fig7n_stats %>%
  mutate(
    x1 = c(1, 2),
    x2 = c(2, 3),
    y_pos = c(
      fig7n_y_range[2] + 0.10 * fig7n_y_span,
      fig7n_y_range[2] + 0.22 * fig7n_y_span
    )
  )

fig7n_colors <- c(
  CHOW = "#BFBFBF",
  MASH = "#D73027",
  MASHantiSPP1 = "#CCA961"
)

fig7n_labels <- c(
  CHOW = "Chow",
  MASH = "MASH",
  MASHantiSPP1 = "MASH + anti-OPN"
)

p7N <- ggplot(fig7n_df, aes(group, mean_score, fill = group)) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.4,
    alpha = 0.75
  ) +
  geom_jitter(
    shape = 21,
    color = "black",
    size = 2.8,
    stroke = 0.4,
    width = 0.10,
    height = 0
  ) +
  geom_segment(
    data = fig7n_stats,
    aes(x = x1, xend = x2, y = y_pos, yend = y_pos),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    data = fig7n_stats,
    aes(x = (x1 + x2) / 2, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 3.5
  ) +
  scale_fill_manual(values = fig7n_colors) +
  scale_x_discrete(labels = fig7n_labels) +
  coord_cartesian(
    ylim = c(
      fig7n_y_range[1] - 0.05 * fig7n_y_span,
      fig7n_y_range[2] + 0.35 * fig7n_y_span
    )
  ) +
  labs(
    x = NULL,
    y = "Fibroblast activation score"
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11),
    legend.position = "none"
  )

print(p7N)

write.csv(
  fig7n_df,
  file.path(OUT_DIR, "Figure7N_FB_activation_score_sample_level.csv"),
  row.names = FALSE
)

write.csv(
  fig7n_stats,
  file.path(OUT_DIR, "Figure7N_FB_activation_score_stats.csv"),
  row.names = FALSE
)

ggsave(
  file.path(OUT_DIR, "Figure7N_FB_activation_score.pdf"),
  p7N,
  width = 4.8,
  height = 4.5,
  device = cairo_pdf
)

ggsave(
  file.path(OUT_DIR, "Figure7N_FB_activation_score.png"),
  p7N,
  width = 4.8,
  height = 4.5,
  units = "in",
  dpi = 600
)

# =============================================================================
# 4. Pseudo-bulk DESeq2 (per-sample sums of fibroblast counts)
#    MASH vs Chow: ~ batch + group; Spp1 KO vs flox (single batch): ~ group
# =============================================================================
meta <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))
counts <- GetAssayData(FB, assay = "RNA", layer = "counts")
pb_counts <- sapply(meta$orig.ident, function(s)
  Matrix::rowSums(counts[, FB$orig.ident == s, drop = FALSE]))
pb_meta <- data.frame(group = meta$Group, batch = meta$Batch, row.names = meta$orig.ident)

run_pb_deseq2 <- function(numerator, denominator, name) {
  keep <- pb_meta$group %in% c(numerator, denominator)
  m  <- pb_counts[, keep]
  md <- pb_meta[keep, ]
  md$group <- factor(md$group, levels = c(denominator, numerator))
  md$batch <- droplevels(factor(md$batch))

  m <- m[filterByExpr(m, group = md$group, min.count = 10, min.total.count = 15), ]
  design <- if (nlevels(md$batch) > 1) ~ batch + group else ~ group

  dds <- DESeq(DESeqDataSetFromMatrix(round(m), md, design), quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("group", numerator, denominator), alpha = 0.05)) %>%
    tibble::rownames_to_column("gene") %>% arrange(padj)

  out <- file.path(OUT_DIR, "pseudobulk_DESeq2", name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(res, file.path(out, "DESeq2_full_results.csv"), row.names = FALSE)
  write.csv(filter(res, !is.na(padj), padj < 0.05), file.path(out, "DESeq2_FDR0.05_DEGs.csv"),
            row.names = FALSE)
  message(name, ": ", sum(res$padj < 0.05, na.rm = TRUE), " genes FDR < 0.05")
  res
}
de_mash <- run_pb_deseq2("MASH", "CHOW", "MASH_vs_CHOW")
de_ko   <- run_pb_deseq2("MASH_SPP1KO", "MASH_flox", "SPP1KO_vs_flox")

# Quadrant plot: genes significant in both contrasts (FDR < 0.05, |log2FC| > 0.25)
fc_cutoff <- 0.25; lim <- 2.5
blacklist <- unique(c(grep("^mt-", rownames(FB), value = TRUE),
                      grep("^(Rpl|Rps|Mrpl|Mrps)", rownames(FB), value = TRUE),
                      grep("^(Hba|Hbb|Alas2)", rownames(FB), value = TRUE),
                      c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty")))
fb_quad <- inner_join(
  de_mash %>% transmute(gene, x = log2FoldChange, padj_x = padj),
  de_ko   %>% transmute(gene, y = log2FoldChange, padj_y = padj), by = "gene") %>%
  filter(!gene %in% blacklist, !is.na(x), !is.na(y)) %>%
  mutate(sig = abs(x) > fc_cutoff & abs(y) > fc_cutoff &
               !is.na(padj_x) & padj_x < 0.05 & !is.na(padj_y) & padj_y < 0.05,
         quadrant = case_when(!sig ~ "NS",
                              x > 0 & y > 0 ~ "Q1: Both Up",
                              x < 0 & y > 0 ~ "Q2: X-Down / Y-Up",
                              x < 0 & y < 0 ~ "Q3: Both Down",
                              TRUE          ~ "Q4: X-Up / Y-Down"),
         dist_val = sqrt(x^2 + y^2))
write.csv(filter(fb_quad, quadrant != "NS") %>% arrange(quadrant, desc(dist_val)),
          file.path(OUT_DIR, "pseudobulk_DESeq2", "quadrant_genes.csv"), row.names = FALSE)

quad_colors <- c("Q1: Both Up" = "#D73027", "Q2: X-Down / Y-Up" = "#4575B4",
                 "Q3: Both Down" = "#1A9850", "Q4: X-Up / Y-Down" = "#F46D43", NS = "grey80")
quad_n <- data.frame(quadrant = names(quad_colors)[1:4],
                     x = c(1, -1, -1, 1) * lim * 0.8, y = c(1, 1, -1, -1) * lim * 0.9) %>%
  mutate(n = sapply(quadrant, function(q) sum(fb_quad$quadrant == q)))
label_genes <- fb_quad %>% filter(quadrant != "NS") %>% group_by(quadrant) %>%
  slice_max(dist_val, n = 8, with_ties = FALSE) %>% ungroup()

p <- ggplot(fb_quad, aes(x, y)) +
  annotate("rect", xmin = c(0, -lim, -lim, 0), xmax = c(lim, 0, 0, lim),
           ymin = c(0, 0, -lim, -lim), ymax = c(lim, lim, 0, 0),
           fill = quad_colors[1:4], alpha = 0.03) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(data = filter(fb_quad, quadrant == "NS"), color = "grey80", size = 0.5, alpha = 0.4) +
  geom_point(data = filter(fb_quad, quadrant != "NS"), aes(color = quadrant), size = 1.8, alpha = 0.8) +
  geom_text_repel(data = label_genes, aes(label = gene, color = quadrant), size = 2.8,
                  max.overlaps = 30, segment.size = 0.3, segment.alpha = 0.5,
                  fontface = "italic", show.legend = FALSE, seed = 123) +
  geom_text(data = quad_n, aes(label = paste0("n=", n), color = quadrant),
            size = 4, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = quad_colors, name = "Quadrant") +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = expression(log[2] * FC ~ "(MASH vs Chow)"),
       y = expression(log[2] * FC ~ "(Spp1 KO vs flox)"),
       title = "Pseudobulk DEGs in fibroblasts") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(OUT_DIR, "pseudobulk_quadplot.pdf"), p, width = 7, height = 7)

saveRDS(FB, file.path(DATA_DIR, "FB_subclusters.rds"))
sessionInfo()
