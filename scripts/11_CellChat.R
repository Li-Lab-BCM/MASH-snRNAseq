# =============================================================================
# CellChat among MP1 (macrophage sub-cluster), FB1 (fibroblast sub-cluster),
# atrial cardiomyocytes and endothelial cells, run per group and compared:
#   MASH vs CHOW, anti-OPN vs MASH, SPP1KO vs flox.
# "Rescue" ligand-receptor pairs: changed in MASH vs CHOW and reversed in
# MASH_SPP1KO vs MASH_flox.
#
# Input : <DATA_DIR>/seurat_final.rds
#         <DATA_DIR>/MPDC_subclusters.rds
#         <DATA_DIR>/FB_subclusters.rds
# Output: <DATA_DIR>/results/CellChat/
# =============================================================================


library(Seurat)
library(CellChat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

DATA_DIR <- "/path/to/MASHproject"
OUT_DIR  <- file.path(DATA_DIR, "results", "CellChat")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(1234)

seurat_object <- readRDS(file.path(DATA_DIR, "seurat_final.rds"))

cols.use <- c("MP1"                  = "#E7298A",
              "FB1"                  = "#FF7F00",
              "Atrial Cardiomyocyte" = "#8491B4FF",
              "Endothelial cell"     = "#3C5488FF")
all_ct <- names(cols.use)

# -----------------------------------------------------------------------------
# 1. Reuse the finalized MP/DC and fibroblast subcluster annotations
# -----------------------------------------------------------------------------
mp <- readRDS(file.path(DATA_DIR, "MPDC_subclusters.rds"))
fb <- readRDS(file.path(DATA_DIR, "FB_subclusters.rds"))

if (!"celltype" %in% colnames(mp@meta.data)) {
  stop("MPDC_subclusters.rds does not contain the finalized MP/DC subtype labels.")
}
if (!"fb_cluster" %in% colnames(fb@meta.data)) {
  stop("FB_subclusters.rds does not contain the finalized fibroblast subtype labels.")
}

mp1_cells <- intersect(
  colnames(mp)[as.character(mp$celltype) == "MP1"],
  colnames(seurat_object)
)
fb1_cells <- intersect(
  colnames(fb)[as.character(fb$fb_cluster) == "FB1"],
  colnames(seurat_object)
)

if (length(mp1_cells) == 0) stop("No MP1 nuclei matched seurat_final.rds.")
if (length(fb1_cells) == 0) stop("No FB1 nuclei matched seurat_final.rds.")

seurat_object$cellchat_group <- NA_character_
seurat_object$cellchat_group[mp1_cells] <- "MP1"
seurat_object$cellchat_group[fb1_cells] <- "FB1"
seurat_object$cellchat_group[
  seurat_object$celltype == "Atrial Cardiomyocyte"
] <- "Atrial Cardiomyocyte"
seurat_object$cellchat_group[
  seurat_object$celltype == "Endothelial cell"
] <- "Endothelial cell"

cat("\nCellChat populations in seurat_final.rds:\n")
print(table(seurat_object$cellchat_group, useNA = "ifany"))

rm(mp, fb, mp1_cells, fb1_cells)

# -----------------------------------------------------------------------------
# 2. CellChat per group
# -----------------------------------------------------------------------------
run_cellchat <- function(group_name) {
  obj <- subset(seurat_object, subset = group == group_name & !is.na(cellchat_group))
  obj <- JoinLayers(obj, assay = "RNA")
  obj$cellchat_group <- factor(obj$cellchat_group, levels = all_ct)

  cc <- createCellChat(object = obj, group.by = "cellchat_group", assay = "RNA")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE,
                          population.size = FALSE, nboot = 100, seed.use = 1)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc)
  liftCellChat(cc, all_ct)
}

groups <- c("CHOW", "MASH", "MASHantiSPP1", "MASH_flox", "MASH_SPP1KO")

missing_groups <- setdiff(groups, unique(as.character(seurat_object$group)))
if (length(missing_groups) > 0) {
  stop("Missing experimental group(s): ", paste(missing_groups, collapse = ", "))
}

group_counts <- table(
  seurat_object$group[!is.na(seurat_object$cellchat_group)],
  seurat_object$cellchat_group[!is.na(seurat_object$cellchat_group)]
)
cat("\nNuclei entering CellChat by group and population:\n")
print(group_counts)

cc_list <- setNames(lapply(groups, run_cellchat), groups)
saveRDS(cc_list, file.path(OUT_DIR, "cellchat_by_group.rds"))

comparisons <- list(
  MASH_vs_CHOW      = c("CHOW", "MASH"),
  AntiOPN_vs_MASH   = c("MASH", "MASHantiSPP1"),
  SPP1KO_vs_flox    = c("MASH_flox", "MASH_SPP1KO")
)

group_display <- c(
  CHOW = "Chow",
  MASH = "MASH",
  MASHantiSPP1 = "anti-OPN",
  MASH_flox = "Spp1 f/f",
  MASH_SPP1KO = "Spp1 f/f;Alb-Cre"
)

# -----------------------------------------------------------------------------
# 3. Comparison plots: interaction number/strength, heatmaps and MP1-to-target bubbles
# -----------------------------------------------------------------------------
for (nm in names(comparisons)) {
  grp <- comparisons[[nm]]
  merged <- mergeCellChat(cc_list[grp], add.names = grp)

  p <- compareInteractions(merged, show.legend = FALSE, group = c(1, 2)) +
    compareInteractions(merged, show.legend = FALSE, group = c(1, 2), measure = "weight")
  ggsave(file.path(OUT_DIR, paste0(nm, "_interactions.pdf")), p, width = 6, height = 4)

  for (m in c("count", "weight")) {
    pdf(file.path(OUT_DIR, paste0(nm, "_heatmap_", m, ".pdf")), width = 3.5, height = 3)
    draw(netVisual_heatmap(merged, measure = m, color.use = cols.use))
    dev.off()
  }

  targets <- c("FB1", "Atrial Cardiomyocyte", "Endothelial cell")
  p_up <- netVisual_bubble(merged, sources.use = "MP1", targets.use = targets,
                           comparison = c(1, 2), max.dataset = 2, angle.x = 45,
                           remove.isolate = TRUE,
                           title.name = paste0("Increased in ", group_display[grp[2]], " (MP1 -> targets)"))
  p_down <- netVisual_bubble(merged, sources.use = "MP1", targets.use = targets,
                             comparison = c(1, 2), max.dataset = 1, angle.x = 45,
                             remove.isolate = TRUE,
                             title.name = paste0("Decreased in ", group_display[grp[2]], " (MP1 -> targets)"))
  ggsave(file.path(OUT_DIR, paste0(nm, "_MP1_bubble.pdf")), p_up + p_down, width = 14, height = 8)
}

# -----------------------------------------------------------------------------
# 4. Rescue ligand-receptor pairs per sender (MASH vs CHOW reversed by SPP1 KO)
#    up rescue  : increased in MASH, decreased in SPP1KO vs flox
#    down rescue: decreased in MASH, increased in SPP1KO vs flox
# -----------------------------------------------------------------------------
lr_cols <- c("source", "target", "ligand", "receptor", "interaction_name", "pathway_name")

get_lr <- function(cc, source_use) {
  subsetCommunication(cc, sources.use = source_use, targets.use = all_ct) %>%
    dplyr::select(all_of(lr_cols), prob, pval)
}

compare_lr <- function(df1, df2, n1, n2, fc_name) {
  full_join(df1 %>% rename(!!paste0("prob_", n1) := prob, !!paste0("pval_", n1) := pval),
            df2 %>% rename(!!paste0("prob_", n2) := prob, !!paste0("pval_", n2) := pval),
            by = lr_cols) %>%
    mutate(across(starts_with("prob_"), ~ replace_na(.x, 0)),
           across(starts_with("pval_"), ~ replace_na(.x, 1)),
           !!fc_name := log2((.data[[paste0("prob_", n2)]] + 1e-6) /
                             (.data[[paste0("prob_", n1)]] + 1e-6)))
}

plot_rescue <- function(df, prefix, title, top_n = 20) {
  if (nrow(df) == 0) return(invisible(NULL))
  top <- head(df, top_n)
  write.csv(top, file.path(OUT_DIR, paste0(prefix, "_top.csv")), row.names = FALSE)

  # Chord diagram (link width = interaction probability in MASH)
  chord <- as.data.frame(top[, c("source", "target", "prob_MASH")])
  sectors <- unique(c(chord$source, chord$target))
  pdf(file.path(OUT_DIR, paste0(prefix, "_chord.pdf")), width = 8, height = 7)
  circos.clear()
  circos.par(gap.after = c(rep(3, length(sectors) - 1), 20), start.degree = 90,
             canvas.xlim = c(-1.55, 1.55), canvas.ylim = c(-1.55, 1.55))
  pos <- chordDiagram(chord, grid.col = cols.use[sectors],
                      col = adjustcolor(cols.use[chord$target], alpha.f = 0.45),
                      directional = 1, direction.type = "arrows",
                      link.arr.type = "big.arrow", annotationTrack = "grid",
                      preAllocateTracks = list(list(track.height = 0.15)),
                      link.sort = TRUE, link.decreasing = TRUE)
  pos$label <- top$interaction_name
  pos$label_x <- pos$x2 - abs(pos$value2) / 2          # label at receptor end
  circos.track(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    d <- pos[pos$cn == CELL_META$sector.index, ]
    if (nrow(d) > 0) circos.text(d$label_x, CELL_META$ylim[2] - 0.05, d$label,
                                 facing = "clockwise", niceFacing = TRUE,
                                 adj = c(0, 0.5), cex = 0.7, col = "grey20")
  })
  title(main = paste0(title, ": top ", nrow(top), " rescue LR"), cex.main = 1, line = -2)
  legend("right", legend = sectors, fill = cols.use[sectors], border = NA, bty = "n")
  circos.clear()
  dev.off()

  # Dot plot of log2 fold changes
  plot_df <- top %>%
    mutate(lr = factor(paste0(interaction_name, " -> ", target),
                       levels = rev(paste0(interaction_name, " -> ", target)))) %>%
    pivot_longer(c(fc_disease, fc_ko), names_to = "comparison", values_to = "logFC") %>%
    mutate(comparison = ifelse(comparison == "fc_disease", "MASH vs CHOW", "SPP1KO vs Flox"))
  p <- ggplot(plot_df, aes(comparison, lr)) +
    geom_point(aes(size = abs(logFC), color = logFC)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "log2FC") +
    scale_size_continuous(range = c(2, 8), name = "|log2FC|") +
    labs(title = paste0(title, ": top ", nrow(top), " rescue LR"), x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 8), panel.grid.major.x = element_blank())
  ggsave(file.path(OUT_DIR, paste0(prefix, "_dotplot.pdf")), p,
         width = 6, height = max(6, nrow(top) * 0.35))
}

for (src in c("MP1", "FB1", "Atrial Cardiomyocyte")) {
  tag <- ifelse(src == "Atrial Cardiomyocyte", "ACM", src)

  cmp_disease <- compare_lr(get_lr(cc_list$CHOW, src), get_lr(cc_list$MASH, src),
                            "CHOW", "MASH", "fc_disease")
  cmp_ko      <- compare_lr(get_lr(cc_list$MASH_flox, src), get_lr(cc_list$MASH_SPP1KO, src),
                            "flox", "KO", "fc_ko")

  rescue_up <- inner_join(filter(cmp_disease, fc_disease > 0) %>% select(all_of(lr_cols), prob_CHOW, prob_MASH, fc_disease),
                          filter(cmp_ko, fc_ko < 0) %>% select(all_of(lr_cols), prob_flox, prob_KO, fc_ko),
                          by = lr_cols) %>%
    mutate(rescue_score = fc_disease - fc_ko) %>% arrange(desc(rescue_score))
  rescue_down <- inner_join(filter(cmp_disease, fc_disease < 0) %>% select(all_of(lr_cols), prob_CHOW, prob_MASH, fc_disease),
                            filter(cmp_ko, fc_ko > 0) %>% select(all_of(lr_cols), prob_flox, prob_KO, fc_ko),
                            by = lr_cols) %>%
    mutate(rescue_score = fc_ko - fc_disease) %>% arrange(desc(rescue_score))

  write.csv(rescue_up,   file.path(OUT_DIR, paste0(tag, "_rescue_up_LR.csv")), row.names = FALSE)
  write.csv(rescue_down, file.path(OUT_DIR, paste0(tag, "_rescue_down_LR.csv")), row.names = FALSE)
  plot_rescue(rescue_up,   paste0(tag, "_rescue_up"),   paste(src, "up rescue"))
  plot_rescue(rescue_down, paste0(tag, "_rescue_down"), paste(src, "down rescue"))
}

sessionInfo()
