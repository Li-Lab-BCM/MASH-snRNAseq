# =============================================================================
# Curated plotting code for manuscript figures.
#
# Included:
#   - barcode-rank plots
#   - Supplemental Figure S7A-D post-QC overview
#   - Figure 3C whole-atrium GSEA
#   - Figure 3E revised MP/DC GSEA
#   - Figure 3L MP/DC Tgfbr1/Cd44
#   - ACM anti-OPN GO:BP plot
#   - ACM Spp1KO Reactome plot
#   - EC anti-OPN GO:BP plot
#   - EC Spp1KO GO:BP plot
#   - Whole-atrium anti-OPN GO:BP plot
#   - MP/DC anti-OPN selected GO:BP plot


library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)


DATA_DIR <- "/path/to/MASHproject"
RAW_DIR <- Sys.getenv("MASH_RAW_DIR", unset = DATA_DIR)
RESULTS_DIR <- file.path(DATA_DIR, "results")
FIG_DIR <- file.path(RESULTS_DIR, "figures")
GSEA_ROOT <- file.path(RESULTS_DIR, "LowExprFiltered_DESeq2_GSEA")

root_dir <- RAW_DIR

sample_names <- c(
  "A1","A2","A3",
  "B1","B2","B3",
  "D1","D2","D3",
  "G1","G2",
  "F1","F3"
)

out_dir <- file.path(FIG_DIR, "BarcodeRank")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


barcode_data <- list()
summary_list <- list()

for (sample in sample_names) {

  cat("\nProcessing:", sample, "\n")

  sample_dir <- file.path(root_dir, sample)

  raw_files <- list.files(
    sample_dir,
    pattern = "(^|_)raw_feature_bc_matrix\\.h5$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(raw_files) == 0) {
    stop("No raw feature-barcode matrix found for: ", sample)
  }

  raw_file <- raw_files[1]

  cat("  ", raw_file, "\n")


  library_name <- basename(raw_file)

  library_name <- sub(
    "_raw_feature_bc_matrix\\.h5$",
    "",
    library_name,
    ignore.case = TRUE
  )


  raw_counts <- Read10X_h5(raw_file)

  if (is.list(raw_counts)) {

    if ("Gene Expression" %in% names(raw_counts)) {
      raw_counts <- raw_counts[["Gene Expression"]]
    } else {
      raw_counts <- raw_counts[[1]]
    }
  }


  umi <- Matrix::colSums(raw_counts)

  umi <- umi[umi > 0]

  umi <- sort(
    umi,
    decreasing = TRUE
  )

  df <- data.frame(
    sample = sample,
    library = library_name,
    rank = seq_along(umi),
    UMI = as.numeric(umi)
  )

  barcode_data[[sample]] <- df

  summary_list[[sample]] <- data.frame(
    sample = sample,
    library = library_name,
    n_raw_barcodes_with_UMI = nrow(df),
    max_UMI = max(df$UMI),
    median_UMI = median(df$UMI)
  )

  rm(raw_counts, umi)
  gc()
}


global_max_rank <- max(
  sapply(barcode_data, function(x) max(x$rank))
)

global_max_umi <- max(
  sapply(barcode_data, function(x) max(x$UMI))
)


x_upper <- 10^ceiling(log10(global_max_rank))
y_upper <- 10^ceiling(log10(global_max_umi))

cat("\nGlobal X maximum:", x_upper, "\n")
cat("Global Y maximum:", y_upper, "\n")


plot_list <- list()

for (sample in sample_names) {

  df <- barcode_data[[sample]]

  library_name <- unique(df$library)

  p <- ggplot(
    df,
    aes(x = rank, y = UMI)
  ) +

    geom_line(
      linewidth = 0.55,
      lineend = "round"
    ) +

    scale_x_log10(
      limits = c(1, x_upper),
      breaks = scales::breaks_log(n = 4),
      labels = scales::label_number(
        big.mark = ","
      ),
      expand = expansion(mult = c(0.02, 0.02))
    ) +

    scale_y_log10(
      limits = c(1, y_upper),
      breaks = scales::breaks_log(n = 5),
      labels = scales::label_number(
        big.mark = ","
      ),
      expand = expansion(mult = c(0.02, 0.02))
    ) +

    labs(
      title = paste0(sample, " (", library_name, ")"),
      x = "Barcode rank",
      y = "Total UMI counts"
    ) +

    theme_classic(
      base_size = 11
    ) +

    theme(
      plot.title = element_text(
        size = 10.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),

      axis.title = element_text(
        size = 10,
        face = "bold"
      ),

      axis.text = element_text(
        size = 8.5
      ),

      axis.line = element_line(
        linewidth = 0.6
      ),

      axis.ticks = element_line(
        linewidth = 0.5
      ),

      plot.margin = margin(
        6, 6, 6, 6
      )
    )

  plot_list[[sample]] <- p
}


combined_plot <- wrap_plots(
  plot_list,
  ncol = 4
) +
  plot_annotation(
    title = "Barcode-rank plots of raw snRNA-seq libraries",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  )


ggsave(
  filename = file.path(
    out_dir,
    "Supplemental_BarcodeRank_13samples_FINAL.pdf"
  ),
  plot = combined_plot,
  width = 14,
  height = 12,
  device = cairo_pdf
)


barcode_summary <- bind_rows(summary_list)

write.csv(
  barcode_summary,
  file.path(
    out_dir,
    "BarcodeRank_summary_13samples.csv"
  ),
  row.names = FALSE
)

print(barcode_summary)

cat(
  "\nDONE\nOutput folder:\n",
  out_dir,
  "\n"
)



# =============================================================================
# Supplemental Figure S7A-D — Post-QC overview
# =============================================================================


library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

qc_rds <- file.path(DATA_DIR, "seurat_final.rds")
if (!file.exists(qc_rds)) {
  stop("Cannot find seurat_final.rds: ", qc_rds)
}

qc_object <- readRDS(qc_rds)
DefaultAssay(qc_object) <- "RNA"

if (!"percent.mt" %in% colnames(qc_object@meta.data)) {
  mt_genes <- grep("^mt-", rownames(qc_object), value = TRUE)
  if (length(mt_genes) == 0) stop("No mitochondrial genes found.")
  qc_object[["percent.mt"]] <- PercentageFeatureSet(qc_object, features = mt_genes)
}

qc_sample_order <- c(
  "CHOW_1", "CHOW_2", "CHOW_3",
  "MASH_1", "MASH_2", "MASH_3",
  "MASHantiSPP1_1", "MASHantiSPP1_2", "MASHantiSPP1_3",
  "MASH_flox_1", "MASH_flox_2",
  "MASH_SPP1KO_1", "MASH_SPP1KO_2"
)

qc_group_map <- c(
  CHOW_1 = "CHOW",
  CHOW_2 = "CHOW",
  CHOW_3 = "CHOW",
  MASH_1 = "MASH",
  MASH_2 = "MASH",
  MASH_3 = "MASH",
  MASHantiSPP1_1 = "MASHantiSPP1",
  MASHantiSPP1_2 = "MASHantiSPP1",
  MASHantiSPP1_3 = "MASHantiSPP1",
  MASH_flox_1 = "MASH_flox",
  MASH_flox_2 = "MASH_flox",
  MASH_SPP1KO_1 = "MASH_SPP1KO",
  MASH_SPP1KO_2 = "MASH_SPP1KO"
)

qc_group_order <- c(
  "CHOW",
  "MASH",
  "MASHantiSPP1",
  "MASH_flox",
  "MASH_SPP1KO"
)

qc_group_cols <- c(
  CHOW = "#BFBFBF",
  MASH = "#FDAFA7",
  MASHantiSPP1 = "#CCA961",
  MASH_flox = "#C95D7B",
  MASH_SPP1KO = "#70C171"
)

qc_group_labels <- c(
  CHOW = "plain('Chow')",
  MASH = "plain('MASH')",
  MASHantiSPP1 = "plain('MASH + anti-OPN')",
  MASH_flox = "italic(Spp1)^{f/f}",
  MASH_SPP1KO = "italic(Spp1)^{f/f}*';'*italic(Alb)^{Cre}"
)

qc_df <- qc_object@meta.data %>%
  mutate(
    nucleus = rownames(.),
    sample = as.character(orig.ident),
    group_qc = unname(qc_group_map[sample])
  ) %>%
  filter(sample %in% qc_sample_order, !is.na(group_qc)) %>%
  mutate(
    sample = factor(sample, levels = qc_sample_order),
    group_qc = factor(group_qc, levels = qc_group_order)
  )

qc_out <- file.path(FIG_DIR, "Supplement_Figure_S7_QC")
dir.create(qc_out, recursive = TRUE, showWarnings = FALSE)

qc_theme <- theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10.5, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "none"
  )

qc_x_scale <- scale_x_discrete(
  breaks = qc_group_order,
  labels = parse(text = unname(qc_group_labels[qc_group_order]))
)

# S7A — nuclei per group after QC
qc_count_group <- qc_df %>%
  count(group_qc, name = "n_nuclei")

p_s7a <- ggplot(qc_count_group, aes(group_qc, n_nuclei, fill = group_qc)) +
  geom_col(width = 0.68, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = comma(n_nuclei)),
    vjust = -0.35, size = 4, color = "black"
  ) +
  scale_fill_manual(values = qc_group_cols) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.15))
  ) +
  qc_x_scale +
  labs(
    title = "Nuclei per group after QC",
    x = NULL,
    y = "Number of nuclei"
  ) +
  qc_theme

plot_qc_group_violin <- function(feature, title_text, ylab) {
  ggplot(qc_df, aes(group_qc, .data[[feature]], fill = group_qc)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      color = "black",
      linewidth = 0.22,
      alpha = 0.88
    ) +
    geom_boxplot(
      width = 0.13,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.3
    ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 23,
      size = 2.6,
      fill = "gold",
      color = "black"
    ) +
    scale_fill_manual(values = qc_group_cols) +
    qc_x_scale +
    labs(
      title = title_text,
      x = NULL,
      y = ylab
    ) +
    qc_theme
}

# S7B-D
p_s7b <- plot_qc_group_violin(
  "nFeature_RNA",
  "Detected genes per nucleus by group",
  "nFeature_RNA"
)

p_s7c <- plot_qc_group_violin(
  "nCount_RNA",
  "UMI counts per nucleus by group",
  "nCount_RNA"
)

p_s7d <- plot_qc_group_violin(
  "percent.mt",
  "Mitochondrial percentage by group",
  "percent.mt (%)"
)

p_s7ad <- (p_s7a | p_s7b) / (p_s7c | p_s7d) +
  plot_annotation(tag_levels = "A")

print(p_s7ad)

ggsave(
  file.path(qc_out, "Supplement_Figure_S7A-D_QC_FINAL.pdf"),
  p_s7ad,
  width = 13.5,
  height = 8.5,
  device = cairo_pdf
)


ggsave(
  file.path(qc_out, "Supplement_Figure_S7A_nuclei_per_group.pdf"),
  p_s7a,
  width = 6.5,
  height = 4.2,
  device = cairo_pdf
)

ggsave(
  file.path(qc_out, "Supplement_Figure_S7B_nFeature_by_group.pdf"),
  p_s7b,
  width = 7.2,
  height = 4.2,
  device = cairo_pdf
)

ggsave(
  file.path(qc_out, "Supplement_Figure_S7C_nCount_by_group.pdf"),
  p_s7c,
  width = 7.2,
  height = 4.2,
  device = cairo_pdf
)

ggsave(
  file.path(qc_out, "Supplement_Figure_S7D_percent_mt_by_group.pdf"),
  p_s7d,
  width = 7.2,
  height = 4.2,
  device = cairo_pdf
)

write.csv(
  qc_count_group,
  file.path(qc_out, "Supplement_Figure_S7A_nuclei_counts.csv"),
  row.names = FALSE
)

rm(qc_object, qc_df)


# 2. Figure 3C — Whole-atrium GSEA
library(dplyr)
library(tibble)
library(ggplot2)

# ============================================================
# Figure 3C
# Whole-atrium preranked GO:BP GSEA
# MASH vs CHOW
# Input: fgseaMultilevel results from script 12
# Representative nonredundant significantly enriched pathways are shown
# Dot color = NES
# Dot size  = -log10(FDR)
# ============================================================

ROOT <- file.path(GSEA_ROOT, "Whole_atrium")
OUT  <- file.path(FIG_DIR, "Figure3C")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

gsea_file <- file.path(
  ROOT,
  "MASH_vs_CHOW",
  "GO_BP_GSEA_lowExprFiltered_DESeq2.csv"
)

if (!file.exists(gsea_file)) stop("File not found: ", gsea_file)

gsea <- read.csv(gsea_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  as_tibble()

if (!"FDR" %in% names(gsea)) {
  if ("padj" %in% names(gsea)) {
    gsea$FDR <- gsea$padj
  } else {
    stop("No FDR or padj column found.")
  }
}


# Diagnostic ranking: strongest enriched pathways by NES
gsea_sig <- gsea %>%
  filter(FDR < 0.05, abs(NES) > 1.2)

top5_up_by_NES <- gsea_sig %>%
  filter(NES > 0) %>%
  arrange(desc(NES), FDR) %>%
  slice_head(n = 5)

top5_down_by_NES <- gsea_sig %>%
  filter(NES < 0) %>%
  arrange(NES, FDR) %>%
  slice_head(n = 5)

cat("\n===== Top 5 UP by NES =====\n")
print(top5_up_by_NES %>% select(pathway, NES, FDR), n = Inf)

cat("\n===== Top 5 DOWN by NES =====\n")
print(top5_down_by_NES %>% select(pathway, NES, FDR), n = Inf)

write.csv(
  bind_rows(
    top5_up_by_NES %>% mutate(Direction = "Up"),
    top5_down_by_NES %>% mutate(Direction = "Down")
  ),
  file.path(OUT, "Figure3C_top5_by_NES_diagnostic.csv"),
  row.names = FALSE
)

# Final display set: representative, nonredundant pathways selected from
# significantly enriched GO:BP terms. These are not labeled as the automatic top 5.
up_paths <- c(
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE",
  "GOBP_INTERLEUKIN_6_PRODUCTION",
  "GOBP_T_CELL_MEDIATED_IMMUNITY",
  "GOBP_CELL_KILLING",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"
)

down_paths <- c(
  "GOBP_CARDIAC_MUSCLE_CELL_CONTRACTION",
  "GOBP_ACTIN_MEDIATED_CELL_CONTRACTION",
  "GOBP_CARDIAC_MUSCLE_CELL_ACTION_POTENTIAL",
  "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY",
  "GOBP_MEMBRANE_REPOLARIZATION"
)

label_map <- c(
  GOBP_ADAPTIVE_IMMUNE_RESPONSE = "Adaptive immune response",
  GOBP_INTERLEUKIN_6_PRODUCTION = "Interleukin-6 production",
  GOBP_T_CELL_MEDIATED_IMMUNITY = "T cell mediated immunity",
  GOBP_CELL_KILLING = "Cell killing",
  GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION = "Antigen processing and presentation",
  GOBP_CARDIAC_MUSCLE_CELL_CONTRACTION = "Cardiac muscle cell contraction",
  GOBP_ACTIN_MEDIATED_CELL_CONTRACTION = "Actin-mediated cell contraction",
  GOBP_CARDIAC_MUSCLE_CELL_ACTION_POTENTIAL = "Cardiac muscle cell action potential",
  GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY = "Positive regulation of synapse assembly",
  GOBP_MEMBRANE_REPOLARIZATION = "Membrane repolarization"
)

selected_paths <- c(up_paths, down_paths)

missing_paths <- setdiff(selected_paths, gsea$pathway)
if (length(missing_paths) > 0) {
  stop(
    "Selected pathway(s) not found in the fgsea result:\n",
    paste(missing_paths, collapse = "\n")
  )
}

plot_df <- gsea %>%
  filter(pathway %in% selected_paths) %>%
  distinct(pathway, .keep_all = TRUE) %>%
  mutate(
    Direction = ifelse(pathway %in% up_paths, "Up", "Down"),
    Selection = "Representative nonredundant pathway",
    label = unname(label_map[pathway]),
    logFDR = -log10(pmax(FDR, .Machine$double.xmin))
  )

direction_fail <- plot_df %>%
  filter((Direction == "Up" & NES <= 0) | (Direction == "Down" & NES >= 0))

threshold_fail <- plot_df %>%
  filter(FDR >= 0.05 | abs(NES) <= 1.2)

if (nrow(direction_fail) > 0) {
  print(direction_fail %>% select(pathway, Direction, NES, FDR))
  stop("One or more selected pathways have an unexpected enrichment direction.")
}

if (nrow(threshold_fail) > 0) {
  print(threshold_fail %>% select(pathway, Direction, NES, FDR))
  stop("One or more selected pathways do not meet FDR < 0.05 and |NES| > 1.2.")
}

display_order <- c(
  "Cardiac muscle cell contraction",
  "Actin-mediated cell contraction",
  "Cardiac muscle cell action potential",
  "Positive regulation of synapse assembly",
  "Membrane repolarization",
  "Adaptive immune response",
  "Interleukin-6 production",
  "T cell mediated immunity",
  "Cell killing",
  "Antigen processing and presentation"
)

plot_df <- plot_df %>%
  mutate(
    Direction = factor(Direction, levels = c("Down", "Up")),
    label = factor(label, levels = rev(display_order))
  )

cat("\n===== Figure 3C selected pathways =====\n")
print(
  plot_df %>%
    arrange(Direction, label) %>%
    select(Direction, Selection, pathway, label, NES, FDR, logFDR),
  n = Inf
)

write.csv(
  plot_df,
  file.path(OUT, "Figure3C_selected_pathways.csv"),
  row.names = FALSE
)

p3C <- ggplot(plot_df, aes(x = Direction, y = label)) +
  geom_point(
    aes(size = logFDR, fill = NES),
    shape = 21,
    color = "black",
    stroke = 0.5
  ) +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size_continuous(
    range = c(3.5, 9),
    name = expression(-log[10](FDR))
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

print(p3C)

ggsave(
  file.path(OUT, "Figure3C_WholeAtrium_GSEA_FINAL.pdf"),
  p3C, width = 9, height = 5.8, device = cairo_pdf
)

# =============================================================================
# 3. Figure 3E revised — MP/DC selected pathway alterations
# =============================================================================


library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggnewscale)

# ============================================================
# Figure 3E
# Refined MP/DC pseudobulk DESeq2-Wald GSEA
#
# Top: MASH vs Chow, GO:BP
# Bottom: Spp1KO vs flox, Reactome
#
# Dot size = -log10(FDR)
# Dot color = NES
# ============================================================

ROOT <- file.path(GSEA_ROOT, "MP_DC")
OUT <- file.path(FIG_DIR, "Figure3E")

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. Read GSEA results
# ============================================================

read_gsea <- function(file) {

  x <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ) %>%
    as_tibble()

  if (!"FDR" %in% names(x)) {
    if ("padj" %in% names(x)) {
      x$FDR <- x$padj
    } else {
      stop("No FDR or padj column found.")
    }
  }

  x %>%
    mutate(
      logFDR = -log10(
        pmax(
          FDR,
          .Machine$double.xmin
        )
      )
    )
}


mash_go <- read_gsea(
  file.path(
    ROOT,
    "MASH_vs_CHOW",
    "GO_BP_GSEA_lowExprFiltered_DESeq2.csv"
  )
)

ko_re <- read_gsea(
  file.path(
    ROOT,
    "SPP1KO_vs_flox",
    "REACTOME_GSEA_lowExprFiltered_DESeq2.csv"
  )
)


# ============================================================
# 2. Final IL-1 pathway used in Figure 3E
# ============================================================

IL1_TERM <- "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_1_PRODUCTION"


# ============================================================
# 3. MASH vs Chow terms
# ============================================================

mash_terms <- c(

  "GOBP_MACROPHAGE_ACTIVATION",

  "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",

  IL1_TERM,

  "GOBP_LEUKOCYTE_CHEMOTAXIS",

  "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN"
)


mash_df <- mash_go %>%
  filter(
    pathway %in% mash_terms,
    NES > 1.2,
    FDR < 0.25
  )


if (nrow(mash_df) != 5) {

  cat("\nMissing MASH pathway(s):\n")

  print(
    setdiff(
      mash_terms,
      mash_df$pathway
    )
  )

  stop(
    "MASH panel does not contain exactly 5 valid pathways."
  )
}


# ============================================================
# 4. Spp1KO vs flox terms
# ============================================================

ko_terms <- c(

  "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",

  "REACTOME_GBP_MEDIATED_HOST_DEFENSE",

  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",

  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",

  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE"
)


ko_df <- ko_re %>%
  filter(
    pathway %in% ko_terms,
    NES < -1.2,
    FDR < 0.25
  )


if (nrow(ko_df) != 5) {

  cat("\nMissing KO pathway(s):\n")

  print(
    setdiff(
      ko_terms,
      ko_df$pathway
    )
  )

  stop(
    "KO panel does not contain exactly 5 valid pathways."
  )
}


# ============================================================
# 5. Labels
# ============================================================

label_map <- c(

  "GOBP_MACROPHAGE_ACTIVATION" =
    "Macrophage activation",

  "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION" =
    "TNF superfamily cytokine production",

  "GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_1_PRODUCTION" =
    "Positive regulation of IL-1 production",

  "GOBP_LEUKOCYTE_CHEMOTAXIS" =
    "Leukocyte chemotaxis",

  "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN" =
    "Response to molecule of bacterial origin",

  "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP" =
    "Senescence-associated secretory phenotype (SASP)",

  "REACTOME_GBP_MEDIATED_HOST_DEFENSE" =
    "GBP-mediated host defense",

  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT" =
    "Aerobic respiration and respiratory electron transport",

  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT" =
    "Respiratory electron transport",

  "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE" =
    "SRP-dependent cotranslational protein targeting"
)


mash_df <- mash_df %>%
  mutate(
    label = unname(
      label_map[pathway]
    ),
    group = "Up in MASH (vs Chow)"
  )


ko_df <- ko_df %>%
  mutate(
    label = unname(
      label_map[pathway]
    ),
    group = "Down in Spp1 KO (vs flox)"
  )


# 6. FORCE EXACT ORDER + NUMERIC Y POSITION
# ============================================================

mash_order <- c(
  "Macrophage activation",
  "TNF superfamily cytokine production",
  unname(label_map[IL1_TERM]),
  "Leukocyte chemotaxis",
  "Response to molecule of bacterial origin"
)

ko_order <- c(
  "Senescence-associated secretory phenotype (SASP)",
  "GBP-mediated host defense",
  "Aerobic respiration and respiratory electron transport",
  "Respiratory electron transport",
  "SRP-dependent cotranslational protein targeting"
)

# final top -> bottom order
final_order <- c(
  mash_order,
  ko_order
)

plot_df <- bind_rows(
  mash_df,
  ko_df
)

# assign y positions:
# top pathway = 10
# bottom pathway = 1
y_lookup <- setNames(
  rev(seq_along(final_order)),
  final_order
)

plot_df <- plot_df %>%
  mutate(
    y_pos = unname(y_lookup[label]),
    x_pos = 1
  )

# safety check
if (any(is.na(plot_df$y_pos))) {
  stop(
    "Some labels were not matched to y positions:\n",
    paste(
      plot_df$label[is.na(plot_df$y_pos)],
      collapse = "\n"
    )
  )
}

cat("\n====================================\n")
cat("FINAL FIGURE 3E\n")
cat("====================================\n")

print(
  plot_df %>%
    arrange(desc(y_pos)) %>%
    select(
      group,
      label,
      NES,
      FDR,
      logFDR,
      y_pos
    ) %>%
    as_tibble(),
  n = Inf
)

write.csv(
  plot_df,
  file.path(
    OUT,
    "Figure3E_selected_pathways.csv"
  ),
  row.names = FALSE
)


up_df <- plot_df %>%
  filter(
    group == "Up in MASH (vs Chow)"
  )

down_df <- plot_df %>%
  filter(
    group == "Down in Spp1 KO (vs flox)"
  )


nes_lim <- max(abs(plot_df$NES), na.rm = TRUE)
nes_lim <- ceiling(nes_lim * 10) / 10

p <- ggplot(plot_df, aes(x = x_pos, y = y_pos)) +

  # central outline
  annotate(
    "rect",
    xmin = 0.84, xmax = 1.16,
    ymin = 0.5, ymax = 10.5,
    fill = NA, color = "grey45", linewidth = 0.55
  ) +

  # divider
  annotate(
    "segment",
    x = 0.84, xend = 1.16,
    y = 5.5, yend = 5.5,
    color = "grey70", linewidth = 0.5
  ) +

  # pathway dots
  geom_point(
    aes(size = logFDR, fill = NES),
    shape = 21,
    color = "grey25",
    stroke = 0.45
  ) +

  # NES: blue = negative, red = positive
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-nes_lim, nes_lim),
    name = "NES"
  ) +

  # significance
  scale_size_continuous(
    range = c(4.5, 10),
    name = expression(-log[10](FDR))
  ) +

  # pathway labels
  scale_y_continuous(
    breaks = rev(seq_along(final_order)),
    labels = final_order,
    limits = c(-0.3, 11.4),
    expand = c(0, 0)
  ) +

  annotate(
    "text",
    x = 1, y = 11.05,
    label = "Up in MASH (vs Chow)",
    fontface = "bold", size = 4.8
  ) +

  annotate(
    "text",
    x = 1, y = -0.05,
    label = "Down in Spp1 KO (vs flox)",
    fontface = "bold", size = 4.8
  ) +

  scale_x_continuous(
    limits = c(0.65, 1.62),
    expand = c(0, 0)
  ) +

  labs(
    title = "Selected pathway alterations in MP/DCs",
    x = NULL,
    y = NULL
  ) +

  coord_cartesian(clip = "off") +

  theme_classic(
    base_size = 11,
    base_family = "Arial"
  ) +

  theme(
    plot.title = element_text(face = "bold", size = 13.5, hjust = 0),
    axis.text.y = element_text(size = 10.5, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(15, 18, 20, 8)
  )

print(p)


ggsave(
  file.path(
    OUT,
    "Figure3E_revised.pdf"
  ),
  p,
  width = 8.7,
  height = 6.8,
  device = cairo_pdf
)

cat(
  "\nSaved to:\n",
  OUT,
  "\n"
)

# =============================================================================
# 4. Figure 3L — MP/DC Tgfbr1 and Cd44
# =============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

rds_file <- file.path(DATA_DIR, "MPDC_subclusters.rds")
out_dir <- file.path(FIG_DIR, "Figure3L_MPDC")
dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
obj <- readRDS(rds_file)
DefaultAssay(obj) <- "RNA"
cat("Reductions:\n")
print(names(obj@reductions))
cat("\nGroups:\n")
print(table(obj$group))

if ("umap" %in% names(obj@reductions)) {

  reduction_use <- "umap"

} else {

  umap_candidates <- grep(
    "umap",
    names(obj@reductions),
    value = TRUE,
    ignore.case = TRUE
  )

  if (length(umap_candidates) == 0) {
    stop("No UMAP reduction found.")
  }

  reduction_use <- umap_candidates[1]
}

cat("\nUsing reduction:", reduction_use, "\n")

keep_groups <- c(
  "CHOW",
  "MASH",
  "MASH_flox",
  "MASH_SPP1KO"
)

cells_keep <- rownames(obj@meta.data)[
  obj$group %in% keep_groups
]

obj4 <- subset(
  obj,
  cells = cells_keep
)

cat("\nCells retained:\n")
print(table(obj4$group))

obj4$group_plot <- case_when(

  obj4$group == "CHOW" ~ "Chow",

  obj4$group == "MASH" ~ "MASH",

  obj4$group == "MASH_flox" ~ "Spp1^{f/f}",

  obj4$group == "MASH_SPP1KO" ~ "Spp1^{f/f}*';'*Alb^{Cre}",

  TRUE ~ as.character(obj4$group)
)

obj4$group_plot <- factor(
  obj4$group_plot,
  levels = c(
    "Chow",
    "MASH",
    "Spp1^{f/f}",
    "Spp1^{f/f}*';'*Alb^{Cre}"
  )
)


genes_to_plot <- c(
  "Tgfbr1",
  "Cd44"
)

missing_genes <- setdiff(
  genes_to_plot,
  rownames(obj4)
)

if (length(missing_genes) > 0) {

  stop(
    "Missing genes: ",
    paste(
      missing_genes,
      collapse = ", "
    )
  )
}

emb <- Embeddings(
  obj4,
  reduction = reduction_use
)[, 1:2, drop = FALSE]

colnames(emb) <- c(
  "UMAP1",
  "UMAP2"
)

emb <- as.data.frame(emb) %>%
  tibble::rownames_to_column("cell_id")



expr <- FetchData(
  obj4,
  vars = genes_to_plot,
  layer = "data"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id")

meta <- obj4@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(
    cell_id,
    group,
    group_plot
  )


plot_df <- emb %>%
  left_join(
    meta,
    by = "cell_id"
  ) %>%
  left_join(
    expr,
    by = "cell_id"
  ) %>%
  pivot_longer(
    cols = all_of(genes_to_plot),
    names_to = "gene",
    values_to = "expression"
  )


print(
  plot_df %>%
    group_by(
      gene,
      group_plot
    ) %>%
    summarise(
      n = n(),
      mean_expr = mean(expression),
      max_expr = max(expression),
      .groups = "drop"
    )
)



make_feature_plot <- function(
    df,
    gene_name
) {

  dat <- df %>%
    filter(
      gene == gene_name
    )

  vmax <- as.numeric(
    quantile(
      dat$expression,
      probs = 0.99,
      na.rm = TRUE
    )
  )

  if (
    is.na(vmax) ||
    vmax <= 0
  ) {

    vmax <- max(
      dat$expression,
      na.rm = TRUE
    )
  }

  if (vmax <= 0) {
    vmax <- 1
  }

  dat <- dat %>%
    mutate(
      expr_plot = pmin(
        expression,
        vmax
      )
    )

  ggplot(
    dat,
    aes(
      x = UMAP1,
      y = UMAP2,
      color = expr_plot
    )
  ) +

    geom_point(
      size = 0.32,
      alpha = 0.9
    ) +

    facet_wrap(
      ~group_plot,
      nrow = 1,
      labeller = label_parsed
    ) +

    scale_color_viridis_c(
      option = "plasma",
      limits = c(
        0,
        vmax
      ),
      oob = squish,
      name = "Expression"
    ) +

    labs(
      title = gene_name,
      x = "UMAP1",
      y = "UMAP2"
    ) +

    coord_fixed() +

    theme_classic(
      base_size = 11,
      base_family = "Arial"
    ) +

    theme(

      plot.title = element_text(
        face = "italic",
        size = 15,
        hjust = 0
      ),

      strip.background = element_blank(),

      strip.text = element_text(
        size = 13,
        color = "black"
      ),

      axis.title = element_text(
        size = 10
      ),

      axis.text = element_text(
        size = 8,
        color = "black"
      ),

      axis.ticks = element_line(
        linewidth = 0.3
      ),

      axis.line = element_line(
        linewidth = 0.4
      ),

      legend.title = element_text(
        size = 9
      ),

      legend.text = element_text(
        size = 8
      ),

      panel.spacing = unit(
        0.5,
        "lines"
      )
    )
}



p1 <- make_feature_plot(
  plot_df,
  "Tgfbr1"
)

p2 <- make_feature_plot(
  plot_df,
  "Cd44"
)

final_plot <- p1 / p2

print(final_plot)


ggsave(
  file.path(
    out_dir,
    "Figure3L_MPDC_Tgfbr1_Cd44.pdf"
  ),
  final_plot,
  width = 12,
  height = 6.5,
  device = cairo_pdf
)

cat(
  "\nSaved to:\n",
  out_dir,
  "\n"
)


# 5. ACM final reciprocal GSEA plots
#    ACM_GOBP_AntiOPN.pdf / ACM_Reactome_SPP1KO.pdf


library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
ROOT <- GSEA_ROOT
OUT <- file.path(FIG_DIR, "ACM_GSEA")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

read_gsea <- function(pop, db, cmp){
  result_file <- if (db == "GO_BP") {
    "GO_BP_GSEA_lowExprFiltered_DESeq2.csv"
  } else if (db == "Reactome") {
    "REACTOME_GSEA_lowExprFiltered_DESeq2.csv"
  } else {
    stop("Unsupported GSEA database: ", db)
  }

  f <- file.path(ROOT, pop, cmp, result_file)
  if (!file.exists(f)) stop("GSEA result not found: ", f)

  x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"FDR" %in% names(x)) {
    if ("padj" %in% names(x)) {
      x$FDR <- x$padj
    } else {
      stop("No FDR or padj column in: ", f)
    }
  }

  x %>%
    mutate(
      label = pathway %>%
        str_remove("^GOBP_|^REACTOME_") %>%
        str_replace_all("_", " ") %>%
        str_to_sentence()
    )
}


strict_pair <- function(pop, db, treatment){

  cmp2 <- if(treatment=="AntiOPN") {
    "AntiOPN_vs_MASH"
  } else {
    "SPP1KO_vs_flox"
  }

  mash <- read_gsea(pop, db, "MASH_vs_CHOW") %>%
    select(pathway, label, MASH_NES=NES, MASH_FDR=FDR)

  treat <- read_gsea(pop, db, cmp2) %>%
    select(pathway, TREAT_NES=NES, TREAT_FDR=FDR)

  inner_join(mash, treat, by="pathway") %>%
    filter(
      MASH_FDR < 0.25,
      abs(MASH_NES) > 1.2,
      TREAT_FDR < 0.25,
      abs(TREAT_NES) > 1.2,
      MASH_NES * TREAT_NES < 0
    ) %>%
    mutate(
      treatment=treatment,
      direction=ifelse(
        MASH_NES > 0,
        "MASH up → treatment down",
        "MASH down → treatment up"
      )
    )
}


plot_pair <- function(x, terms, title, treat_label){

  d <- x %>%
    filter(label %in% terms) %>%
    mutate(order=match(label, terms)) %>%
    arrange(order)

  if(nrow(d)==0) stop("None of the selected pathways were found.")

  long <- bind_rows(
    d %>%
      transmute(
        label, order,
        comparison="MASH vs Chow",
        NES=MASH_NES,
        FDR=MASH_FDR
      ),
    d %>%
      transmute(
        label, order,
        comparison=treat_label,
        NES=TREAT_NES,
        FDR=TREAT_FDR
      )
  ) %>%
    mutate(
      comparison=factor(
        comparison,
        levels=c("MASH vs Chow", treat_label)
      ),
      label=factor(
        label,
        levels=rev(unique(d$label))
      ),
      logFDR=-log10(pmax(FDR, 1e-300))
    )

  ggplot(long, aes(comparison, label)) +
    geom_point(
      aes(size=logFDR, fill=NES),
      shape=21,
      color="black",
      stroke=0.35
    ) +
    scale_fill_gradient2(
      low="#2166AC",
      mid="white",
      high="#B2182B",
      midpoint=0,
      name="NES"
    ) +
    scale_size_continuous(
      range=c(4,9),
      name=expression(-log[10](FDR))
    ) +
    labs(
      x=NULL,
      y=NULL,
      title=title
    ) +
    theme_classic(base_size=12) +
    theme(
      plot.title=element_text(face="bold", hjust=0.5),
      axis.text.x=element_text(angle=25, hjust=1),
      axis.text.y=element_text(size=10),
      legend.position="right"
    )
}

savep <- function(p, name, w=8, h=5.5){
  ggsave(
    file.path(OUT, paste0(name, ".pdf")),
    p, width=w, height=h
  )
  ggsave(
    file.path(OUT, paste0(name, ".png")),
    p, width=w, height=h, dpi=900
  )
}

# ============================================================
# 1. ACM — anti-OPN — GO:BP
# ============================================================
acm_anti <- strict_pair(
  "ACM",
  "GO_BP",
  "AntiOPN"
)

acm_anti_terms <- c(
  "Erk1 and erk2 cascade",
  "Cytokine production",
  "Chemokine production",
  "Unsaturated fatty acid metabolic process",
  "Response to nutrient",
  "Apoptotic signaling pathway"
)

p_acm_anti <- plot_pair(
  acm_anti,
  acm_anti_terms,
  "Reciprocal GO:BP enrichment in ACM",
  "anti-OPN vs MASH"
)

savep(
  p_acm_anti,
  "ACM_GOBP_AntiOPN",
  8, 5
)

# ============================================================
# 2. ACM — Spp1KO — Reactome
# ============================================================
acm_ko <- strict_pair(
  "ACM",
  "Reactome",
  "SPP1KO"
)

acm_ko_terms <- c(
  "Cytokine signaling in immune system",
  "Signaling by interleukins",
  "Interferon signaling",
  "Potassium channels",
  "Inwardly rectifying k channels"
)

p_acm_ko <- plot_pair(
  acm_ko,
  acm_ko_terms,
  "Reciprocal Reactome enrichment in ACM",
  "Spp1 KO vs flox"
)

savep(
  p_acm_ko,
  "ACM_Reactome_SPP1KO",
  8, 5
)

# ============================================================

# =============================================================================
# 6. EC final GO:BP rescue plots
#    EC_AntiOPN_GOBP.pdf / EC_Spp1KO_GOBP.pdf
# =============================================================================


library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)
ROOT <- file.path(GSEA_ROOT, "EC")
MASH_FILE <- file.path(ROOT, "MASH_vs_CHOW", "GO_BP_GSEA_lowExprFiltered_DESeq2.csv")
ANTI_FILE <- file.path(ROOT, "AntiOPN_vs_MASH", "GO_BP_GSEA_lowExprFiltered_DESeq2.csv")
KO_FILE   <- file.path(ROOT, "SPP1KO_vs_flox", "GO_BP_GSEA_lowExprFiltered_DESeq2.csv")
OUT <- file.path(FIG_DIR, "EC_GSEA")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

read_gsea <- function(file){
  x <- read.csv(file,stringsAsFactors=FALSE,check.names=FALSE)

  if(!"FDR" %in% names(x)){
    if("padj" %in% names(x)) x$FDR <- x$padj
    else if("p.adjust" %in% names(x)) x$FDR <- x$p.adjust
  }

  x %>%
    mutate(
      label=pathway %>%
        str_remove("^GOBP_") %>%
        str_replace_all("_"," ") %>%
        str_to_sentence(),
      key=label %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]",""),
      logFDR=-log10(pmax(FDR,1e-300))
    )
}

mash <- read_gsea(MASH_FILE)
anti <- read_gsea(ANTI_FILE)
ko   <- read_gsea(KO_FILE)


shared <- c(
  "Blood vessel morphogenesis",
  "Response to growth factor",
  "Actin filament organization",
  "Cell junction organization",
  "Cell substrate adhesion"
)

anti_emphasized <- c(
  "Focal adhesion assembly",
  "Cell matrix adhesion",
  "Lamellipodium organization"
)

ko_emphasized <- c(
  "Small gtpase mediated signal transduction",
  "Rho protein signal transduction",
  "Cell surface receptor protein tyrosine kinase signaling pathway"
)

anti_paths <- c(shared,anti_emphasized)
ko_paths   <- c(shared,ko_emphasized)

make_key <- function(x){
  str_replace_all(str_to_lower(x),"[^a-z0-9]","")
}


extract_set <- function(mash_df,treat_df,paths,emphasized,treat_name){

  sel <- tibble(
    label_display=paths,
    key=make_key(paths),
    order=seq_along(paths),
    section=ifelse(paths %in% shared,
                   "Shared EC rescue",
                   emphasized)
  )

  mash_x <- mash_df %>%
    select(
      key,
      MASH_NES=NES,
      MASH_FDR=FDR,
      MASH_logFDR=logFDR
    )

  treat_x <- treat_df %>%
    select(
      key,
      TREAT_NES=NES,
      TREAT_FDR=FDR,
      TREAT_logFDR=logFDR
    )

  out <- sel %>%
    left_join(mash_x,by="key") %>%
    left_join(treat_x,by="key")

  cat("\n============================================\n")
  cat(treat_name,"\n")
  cat("============================================\n")

  print(
    out %>%
      select(
        label_display,
        section,
        MASH_NES,
        MASH_FDR,
        TREAT_NES,
        TREAT_FDR
      )
  )

  if(any(is.na(out$MASH_NES) | is.na(out$TREAT_NES))){
    cat("\nNOT FOUND:\n")
    print(
      out %>%
        filter(is.na(MASH_NES) | is.na(TREAT_NES)) %>%
        pull(label_display)
    )
  }

  out
}

anti_sel <- extract_set(
  mash,
  anti,
  anti_paths,
  "Anti-OPN emphasized",
  "EC: Anti-OPN"
)

ko_sel <- extract_set(
  mash,
  ko,
  ko_paths,
  "Spp1KO emphasized",
  "EC: Spp1KO"
)


to_long <- function(x,treat_label,paths){

  x %>%
    pivot_longer(
      cols=c(MASH_NES,TREAT_NES),
      names_to="Comparison",
      values_to="NES"
    ) %>%
    mutate(
      FDR=ifelse(
        Comparison=="MASH_NES",
        MASH_FDR,
        TREAT_FDR
      ),
      logFDR=ifelse(
        Comparison=="MASH_NES",
        MASH_logFDR,
        TREAT_logFDR
      ),
      Comparison=recode(
        Comparison,
        "MASH_NES"="MASH vs CHOW",
        "TREAT_NES"=treat_label
      ),
      Comparison=factor(
        Comparison,
        levels=c("MASH vs CHOW",treat_label)
      ),
      label_display=factor(
        label_display,
        levels=rev(paths)
      )
    )
}

anti_long <- to_long(
  anti_sel,
  "Anti-OPN vs MASH",
  anti_paths
)

ko_long <- to_long(
  ko_sel,
  "Spp1KO vs flox",
  ko_paths
)


make_plot <- function(df,title){

  ggplot(df,aes(Comparison,label_display)) +

    geom_hline(
      yintercept=3.5,
      linewidth=.5,
      linetype="dashed"
    ) +

    geom_point(
      aes(size=logFDR,fill=NES),
      shape=21,
      color="black",
      stroke=.35
    ) +

    scale_fill_gradient2(
      low="#2166AC",
      mid="white",
      high="#B2182B",
      midpoint=0,
      limits=c(-2.2,2.2),
      name="NES"
    ) +

    scale_size_continuous(
      range=c(3,8),
      name=expression(-log[10](FDR))
    ) +

    labs(
      x=NULL,
      y=NULL,
      title=title
    ) +

    theme_classic(base_size=13) +

    theme(
      plot.title=element_text(
        face="bold",
        hjust=.5,
        size=15
      ),
      axis.text.x=element_text(
        face="bold",
        size=11
      ),
      axis.text.y=element_text(
        size=10.5
      ),
      axis.ticks=element_blank(),
      legend.position="right"
    )
}

p_anti <- make_plot(
  anti_long,
  "EC rescue by anti-OPN"
)

p_ko <- make_plot(
  ko_long,
  "EC rescue by Spp1 deletion"
)

p_final <- p_anti + p_ko +
  plot_layout(
    guides="collect",
    widths=c(1,1)
  ) &
  theme(
    legend.position="right"
  )

p_final

ggsave(
  file.path(OUT, "EC_AntiOPN_GOBP.pdf"),
  p_anti, width = 7, height = 6
)

ggsave(
  file.path(OUT, "EC_Spp1KO_GOBP.pdf"),
  p_ko, width = 7, height = 6
)

write.csv(
  anti_sel,
  file.path(OUT, "EC_AntiOPN_selected_pathways.csv"),
  row.names = FALSE
)

write.csv(
  ko_sel,
  file.path(OUT, "EC_Spp1KO_selected_pathways.csv"),
  row.names = FALSE
)

p_anti
p_ko

# =============================================================================
# 7. Whole-atrium anti-OPN GO:BP plot
# =============================================================================


library(dplyr)
library(stringr)
library(ggplot2)
library(tibble)
WHOLE_FILE <- file.path(
  GSEA_ROOT, "Whole_atrium", "AntiOPN_vs_MASH",
  "GO_BP_GSEA_lowExprFiltered_DESeq2.csv"
)
OUT <- file.path(FIG_DIR, "AntiOPN_GSEA")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
read_gsea <- function(file){

  x <- read.csv(file, stringsAsFactors=FALSE, check.names=FALSE)

  if(!"FDR" %in% names(x)){
    if("padj" %in% names(x)) x$FDR <- x$padj
    else if("p.adjust" %in% names(x)) x$FDR <- x$p.adjust
    else stop("No FDR/padj/p.adjust column in: ", file)
  }

  if(!"pathway" %in% names(x)){
    if("Description" %in% names(x)) x$pathway <- x$Description
    else stop("No pathway/Description column in: ", file)
  }

  x %>%
    mutate(
      NES=as.numeric(NES),
      FDR=as.numeric(FDR),
      logFDR=-log10(pmax(FDR,1e-300)),
      label=pathway %>%
        str_remove("^GOBP_") %>%
        str_replace_all("_"," ") %>%
        str_squish() %>%
        str_to_sentence(),
      label_norm=str_to_lower(label)
    )
}

whole_all <- read_gsea(WHOLE_FILE)

# ============================================================
# FIGURE 1
# WHOLE ATRIUM
# Anti-OPN vs MASH
# Top 5 by NES
# ============================================================
whole <- whole_all %>%
  filter(FDR<0.25, abs(NES)>1.2) %>%
  filter(
    !str_detect(
      label_norm,
      "cilium|ciliary|flagell|microtubule"
    )
  )

whole_up <- whole %>%
  filter(NES>0) %>%
  arrange(desc(NES)) %>%
  slice_head(n=5)

whole_down <- whole %>%
  filter(NES<0) %>%
  arrange(NES) %>%
  slice_head(n=5)

cat("\n========================================\n")
cat("WHOLE ATRIUM: TOP 5 UP BY NES\n")
cat("========================================\n")
print(whole_up %>% select(label,NES,FDR))

cat("\n========================================\n")
cat("WHOLE ATRIUM: TOP 5 DOWN BY NES\n")
cat("========================================\n")
print(whole_down %>% select(label,NES,FDR))

whole_plot <- bind_rows(
  whole_up   %>% mutate(Direction="Up"),
  whole_down %>% mutate(Direction="Down")
) %>%
  mutate(
    Direction=factor(Direction,levels=c("Down","Up")),
    label=factor(label,levels=rev(c(whole_up$label,whole_down$label)))
  )

p_whole <- ggplot(whole_plot,aes(Direction,label)) +
  geom_point(
    aes(size=logFDR,fill=NES),
    shape=21,color="black",stroke=.3
  ) +
  scale_fill_gradient2(
    low="#2166AC",mid="white",high="#D73027",
    midpoint=0,name="NES"
  ) +
  scale_size_continuous(
    range=c(3,8),
    name=expression(-log[10](FDR))
  ) +
  labs(
    x=NULL,y=NULL,
    title="Top altered GO:BP pathways",
    subtitle="Anti-OPN vs MASH"
  ) +
  theme_classic(base_size=13) +
  theme(
    plot.title=element_text(face="bold",hjust=.5,size=15),
    plot.subtitle=element_text(hjust=.5,size=12),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=11),
    axis.ticks=element_blank(),
    legend.position="right"
  )

# ============================================================

ggsave(
  file.path(OUT, "WholeAtrium_AntiOPN_Top5_by_NES.pdf"),
  p_whole, width = 8, height = 5.5
)
write.csv(
  whole_plot,
  file.path(OUT, "WholeAtrium_AntiOPN_selected.csv"),
  row.names = FALSE
)

p_whole

# =============================================================================
# 8. MP/DC anti-OPN selected GO:BP plot
# =============================================================================

library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(forcats)
library(scales)
gsea_file <- file.path(
  GSEA_ROOT, "MP_DC", "AntiOPN_vs_MASH", "GO_BP_FDR025_NES12.csv"
)
out_dir <- file.path(FIG_DIR, "MPDC_AntiOPN_GOBP")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
gsea <- read.csv(gsea_file, check.names = FALSE, stringsAsFactors = FALSE)
cat("Columns in GSEA file:\n")
print(colnames(gsea))

find_first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

col_pathway <- find_first_existing(gsea, c("pathway", "ID", "gs_name", "gene_set"))
col_label   <- find_first_existing(gsea, c("label", "Description", "description", "term"))
col_nes     <- find_first_existing(gsea, c("NES", "nes"))
col_fdr     <- find_first_existing(gsea, c("FDR", "padj", "p.adjust", "qvalue", "qvalues"))

if (is.na(col_nes) | is.na(col_fdr)) stop("Cannot find NES or FDR column.")
if (is.na(col_label) & is.na(col_pathway)) stop("Cannot find pathway/label column.")

if (is.na(col_label)) col_label <- col_pathway
if (is.na(col_pathway)) col_pathway <- col_label

plot_df <- gsea %>%
  mutate(
    pathway_raw = .data[[col_pathway]],
    label_raw = .data[[col_label]],
    NES = as.numeric(.data[[col_nes]]),
    FDR = as.numeric(.data[[col_fdr]])
  ) %>%
  mutate(
    label = label_raw %>%
      str_replace("^GOBP[_:]", "") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    minus_log10_FDR = -log10(FDR)
  ) %>%
  filter(!is.na(NES), !is.na(FDR), !is.na(label))

cat("\nPreview:\n")
print(head(plot_df[, c("pathway_raw", "label", "NES", "FDR")]))


down_targets <- c(
  "Positive regulation of non canonical nf kappab signal transduction",
  "Chemokine production",
  "Tumor necrosis factor superfamily cytokine production",
  "Response to endoplasmic reticulum stress",
  "Regulation of proteolysis involved in protein catabolic process"
)

up_targets <- c(
  "Glycerophospholipid metabolic process",
  "Phospholipid metabolic process",
  "Glycerolipid biosynthetic process"
)

plot_df <- plot_df %>% mutate(label_key = str_to_lower(label))
down_key <- str_to_lower(down_targets)
up_key <- str_to_lower(up_targets)

down_df <- plot_df %>%
  filter(label_key %in% down_key) %>%
  distinct(label_key, .keep_all = TRUE) %>%
  mutate(direction = "Down in anti-OPN")

up_df <- plot_df %>%
  filter(label_key %in% up_key) %>%
  distinct(label_key, .keep_all = TRUE) %>%
  mutate(direction = "Up in anti-OPN")

cat("\nMatched DOWN pathways:\n")
print(down_df[, c("label", "NES", "FDR")], row.names = FALSE)

cat("\nMatched UP pathways:\n")
print(up_df[, c("label", "NES", "FDR")], row.names = FALSE)

not_found_down <- down_targets[!down_key %in% down_df$label_key]
not_found_up <- up_targets[!up_key %in% up_df$label_key]

if (length(not_found_down) > 0) {
  cat("\nNOT FOUND in DOWN targets:\n")
  print(not_found_down)
}
if (length(not_found_up) > 0) {
  cat("\nNOT FOUND in UP targets:\n")
  print(not_found_up)
}

down_df <- down_df %>% arrange(NES)
up_df <- up_df %>% arrange(NES)

final_df <- bind_rows(down_df, up_df) %>%
  mutate(
    label = factor(label, levels = rev(c(down_df$label, up_df$label))),
    direction = factor(direction, levels = c("Down in anti-OPN", "Up in anti-OPN"))
  )

p <- ggplot(final_df, aes(x = 1, y = label)) +
  geom_point(aes(size = minus_log10_FDR, color = NES)) +
  scale_color_gradient2(low = "#3B6FB6", mid = "white", high = "#C23B3B", midpoint = 0, name = "NES") +
  scale_size(range = c(3, 10), name = expression(-log[10](FDR))) +
  facet_grid(direction ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "MP/DC: Anti-OPN vs MASH", subtitle = "Selected GO Biological Process pathways") +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.y = element_text(color = "black", size = 11),
    legend.position = "right",
    panel.spacing = unit(1.0, "lines")
  )

print(p)

pdf_file <- file.path(out_dir, "MPDC_AntiOPN_vs_MASH_selected_GOBP_dotplot.pdf")
csv_file <- file.path(out_dir, "MPDC_AntiOPN_vs_MASH_selected_GOBP_dotplot_data.csv")

ggsave(pdf_file, p, width = 12.2, height = 5.8, device = cairo_pdf)
write.csv(final_df[, c("direction", "label", "NES", "FDR", "minus_log10_FDR", "pathway_raw")], csv_file, row.names = FALSE)

cat("\nSaved files:\n")
print(c(pdf_file, csv_file))
sessionInfo()