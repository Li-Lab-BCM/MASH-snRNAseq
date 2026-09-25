# ============================================================
# 12_pseudobulk_DESeq2_GSEA_all_populations.R
#
# Low-expression-filtered pseudobulk DESeq2 + GSEA
# Populations:
#   Whole_atrium  (all nuclei)
#   FB            (Fibroblast)
#   MP_DC         (MP/DC)
#   EC            (Endothelial cell)
#   ACM           (Atrial Cardiomyocyte)
#
# Workflow (per population):
# nuclei of the population
# -> sample-level pseudobulk
# -> low-expression filtering (edgeR::filterByExpr)
# -> DESeq2 full model: ~ batch + group
# -> Wald statistic ranking
# -> preranked GSEA
# -> GO:BP + Reactome
#
# Comparisons:
# 1. MASH vs CHOW
# 2. Anti-OPN vs MASH
# 3. SPP1KO vs flox
#
# INPUT : <DATA_DIR>/seurat_final.rds
# OUTPUT: <DATA_DIR>/results/LowExprFiltered_DESeq2_GSEA/<population>/
# ============================================================


library(Seurat)
library(SeuratObject)
library(Matrix)
library(edgeR)
library(DESeq2)
library(dplyr)
library(tibble)
library(fgsea)
library(msigdbr)

if (requireNamespace("future", quietly = TRUE)) {
  future::plan("sequential")
}

set.seed(1234)
options(width = 250)


# ============================================================
# 0. PATH
# ============================================================

# Set the project data directory before running
DATA_DIR <- "/path/to/MASHproject"
REPO_DIR  <- file.path(DATA_DIR, "MASH-snRNAseq")
INPUT_RDS <- file.path(DATA_DIR, "seurat_final.rds")
OUT_ROOT  <- file.path(DATA_DIR, "results", "LowExprFiltered_DESeq2_GSEA")

dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(INPUT_RDS))


# ============================================================
# 1. POPULATIONS
# NULL = all nuclei (whole atrium)
# ============================================================

populations <- list(
  Whole_atrium = NULL,
  FB           = "Fibroblast",
  MP_DC        = "MP/DC",
  EC           = "Endothelial cell",
  ACM          = "Atrial Cardiomyocyte"
)


# ============================================================
# 2. SAMPLE DEFINITION + SAMPLE METADATA
# From metadata/sample_metadata.tsv (13 libraries)
# ============================================================

sample_tsv <- read.delim(file.path(REPO_DIR, "metadata", "sample_metadata.tsv"))

expected_samples <- sample_tsv$orig.ident

group_label <- c(
  CHOW         = "CHOW",
  MASH         = "MASH",
  MASHantiSPP1 = "AntiOPN",
  MASH_flox    = "flox",
  MASH_SPP1KO  = "SPP1KO"
)

meta <- data.frame(
  sample = expected_samples,
  group  = unname(group_label[sample_tsv$Group]),
  batch  = sample_tsv$Batch,
  row.names = expected_samples,
  stringsAsFactors = FALSE
)

meta$group <- factor(meta$group, levels = c("CHOW", "MASH", "AntiOPN", "flox", "SPP1KO"))
meta$batch <- factor(meta$batch, levels = c("batch1", "batch2", "batch3"))

if (length(expected_samples) != 13) stop("Expected 13 sequencing libraries.")
if (anyDuplicated(expected_samples)) stop("Duplicated sample IDs.")
if (anyNA(meta$group) || anyNA(meta$batch)) stop("Unrecognized group or batch label.")

cat("\n========================================\n")
cat("SAMPLE METADATA\n")
cat("========================================\n")
print(meta)
cat("\nGroup x batch:\n")
print(table(meta$group, meta$batch))


# ============================================================
# 3. CHECK DESIGN MATRIX
# ============================================================

design_mat <- model.matrix(~ batch + group, data = meta)

cat("\nDesign rank:", qr(design_mat)$rank, "/", ncol(design_mat), "\n")

if (qr(design_mat)$rank < ncol(design_mat)) {
  stop("Design matrix is not full rank.")
}


# ============================================================
# 4. THREE COMPARISONS
# ============================================================

contrasts <- list(
  MASH_vs_CHOW    = c("group", "MASH",    "CHOW"),
  AntiOPN_vs_MASH = c("group", "AntiOPN", "MASH"),
  SPP1KO_vs_flox  = c("group", "SPP1KO",  "flox")
)


# ============================================================
# 5. TECHNICAL GENES
# Removed ONLY from GSEA ranking
# ============================================================

sex_remove <- c("Xist", "Tsix", "Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Zfy1", "Zfy2", "Sry")

is_technical <- function(g) {
  grepl("^mt-", g, ignore.case = TRUE) |
    grepl("^Rpl[0-9]", g) |
    grepl("^Rps[0-9]", g) |
    grepl("^Hba[-0-9]", g, ignore.case = TRUE) |
    grepl("^Hbb[-0-9]", g, ignore.case = TRUE) |
    g %in% sex_remove
}


# ============================================================
# 6. MSigDB GO:BP + REACTOME
# MSigDB human collections mapped to mouse orthologs
# ============================================================

get_msig <- function(collection_name, subcollection_name) {
  fm <- names(formals(msigdbr::msigdbr))
  if ("db_species" %in% fm) {                      # msigdbr >= 10
    msigdbr::msigdbr(db_species = "HS", species = "Mus musculus",
                     collection = collection_name, subcollection = subcollection_name)
  } else {                                         # msigdbr < 10: human sets mapped to mouse
    msigdbr::msigdbr(species = "Mus musculus",
                     category = collection_name, subcategory = subcollection_name)
  }
}

get_symbol_col <- function(x) {
  hit <- intersect(c("gene_symbol", "db_gene_symbol", "human_gene_symbol"), colnames(x))
  if (length(hit) == 0) stop("Cannot find symbol column.")
  hit[1]
}

cat("\nLoading GO:BP...\n")
msig_go <- get_msig("C5", "GO:BP")
cat("\nLoading Reactome...\n")
msig_re <- get_msig("C2", "CP:REACTOME")

go_pathways       <- lapply(split(msig_go[[get_symbol_col(msig_go)]], msig_go$gs_name), unique)
reactome_pathways <- lapply(split(msig_re[[get_symbol_col(msig_re)]], msig_re$gs_name), unique)

cat("GO:BP pathways:", length(go_pathways), "\n")
cat("Reactome pathways:", length(reactome_pathways), "\n")


# ============================================================
# 7. fgseaMultilevel
# ============================================================

gsea_min_size <- c(
  Whole_atrium = 10,
  FB           = 15,
  MP_DC        = 10,
  EC           = 15,
  ACM          = 15
)

run_gsea <- function(ranks, pathways, comparison, database, min_size) {

  x <- fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = min_size,
    maxSize = 500,
    eps = 0,
    scoreType = "std",
    nproc = 1
  )

  x <- as.data.frame(x)
  x$leadingEdge <- vapply(x$leadingEdge, function(z) paste(z, collapse = ";"), character(1))

  x %>%
    as_tibble() %>%
    mutate(
      comparison = comparison,
      database = database,
      FDR = padj,
      minus_log10_FDR = -log10(pmax(padj, .Machine$double.xmin))
    ) %>%
    arrange(FDR, desc(abs(NES)))
}


# ============================================================
# 8. LOAD OBJECT + RNA COUNTS
# ============================================================

cat("\nLoading:\n", INPUT_RDS, "\n")
obj <- readRDS(INPUT_RDS)
DefaultAssay(obj) <- "RNA"

counts   <- SeuratObject::LayerData(obj[["RNA"]], layer = "counts")
obj_meta <- obj@meta.data[colnames(counts), c("orig.ident", "celltype"), drop = FALSE]
stopifnot(identical(rownames(obj_meta), colnames(counts)))
stopifnot(setequal(unique(obj_meta$orig.ident), expected_samples))

cat("\nCell-type labels:\n")
print(table(obj_meta$celltype))

rm(obj)
gc()


# ============================================================
# 9. RUN ONE POPULATION
# ============================================================

run_population <- function(pop_name, celltype_label) {

  cat("\n\n############################################################\n")
  cat("POPULATION:", pop_name, "\n")
  cat("############################################################\n")

  out_dir <- file.path(OUT_ROOT, pop_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------
  # 9.1 Nuclei of this population
  # ----------------------------------------------------------
  cells <- if (is.null(celltype_label)) rownames(obj_meta) else
    rownames(obj_meta)[obj_meta$celltype == celltype_label]

  sample_factor <- factor(obj_meta[cells, "orig.ident"], levels = expected_samples)

  cell_number <- data.frame(sample = expected_samples,
                            n_nuclei = as.vector(table(sample_factor)))
  cat("\nNUCLEI PER SAMPLE\n")
  print(cell_number, row.names = FALSE)
  write.csv(cell_number, file.path(out_dir, "nuclei_per_sample.csv"), row.names = FALSE)

  # ----------------------------------------------------------
  # 9.2 Sample-level pseudobulk sum
  # genes x nuclei  %*%  nuclei x samples
  # ----------------------------------------------------------
  sample_design <- Matrix::sparse.model.matrix(~ 0 + sample_factor)
  colnames(sample_design) <- expected_samples

  pb_counts <- as.matrix(counts[, cells] %*% sample_design)
  colnames(pb_counts) <- expected_samples

  lib_sizes <- colSums(pb_counts)
  if (any(lib_sizes == 0)) {
    stop("Zero-count pseudobulk sample(s): ", paste(names(lib_sizes)[lib_sizes == 0], collapse = ", "))
  }

  if (any(abs(pb_counts - round(pb_counts)) > 1e-6)) {
    cat("\nNOTE: Non-integer values detected; rounding pseudobulk counts for DESeq2.\n")
  }
  pb_counts <- round(pb_counts)
  storage.mode(pb_counts) <- "integer"

  write.csv(data.frame(sample = colnames(pb_counts), library_size = colSums(pb_counts)),
            file.path(out_dir, "pseudobulk_library_sizes.csv"), row.names = FALSE)

  # ----------------------------------------------------------
  # 9.3 Low-expression filtering
  # edgeR::filterByExpr with group
  # ----------------------------------------------------------
  keep <- edgeR::filterByExpr(
    pb_counts,
    group = meta$group,
    min.count = 10,
    min.total.count = 15
  )

  cat("\nLOW-EXPRESSION FILTERING\n")
  cat("Genes before filtering:", nrow(pb_counts), "\n")
  cat("Genes retained:", sum(keep), "(", round(100 * mean(keep), 2), "% )\n")

  pb_filt <- pb_counts[keep, , drop = FALSE]
  write.csv(data.frame(gene = rownames(pb_counts), retained = keep),
            file.path(out_dir, "low_expression_filter.csv"), row.names = FALSE)

  # ----------------------------------------------------------
  # 9.4 DESeq2 full model
  # ----------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = pb_filt, colData = meta, design = ~ batch + group)
  dds <- DESeq(dds, quiet = TRUE)
  saveRDS(dds, file.path(out_dir, "DESeq2_fullmodel_dds.rds"))

  # ----------------------------------------------------------
  # 9.5 DESeq2 results + Wald ranks + GSEA
  # ----------------------------------------------------------
  gsea_results <- list()
  deg_counts   <- list()

  for (cmp in names(contrasts)) {

    cat("\n==========", pop_name, "|", cmp, "==========\n")

    cmp_dir <- file.path(out_dir, cmp)
    dir.create(cmp_dir, recursive = TRUE, showWarnings = FALSE)

    res_df <- as.data.frame(results(dds, contrast = contrasts[[cmp]], alpha = 0.05)) %>%
      rownames_to_column("gene") %>%
      as_tibble()
    write.csv(res_df, file.path(cmp_dir, "DESeq2_results.csv"), row.names = FALSE)

    n_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
    n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
    cat("DEGs FDR < 0.05:", n_up + n_down, "| Up:", n_up, "| Down:", n_down, "\n")
    deg_counts[[cmp]] <- c(up = n_up, down = n_down)

    # GSEA ranking = DESeq2 Wald statistic
    rank_df <- res_df %>%
      filter(!is.na(stat), is.finite(stat), !is_technical(gene)) %>%
      distinct(gene, .keep_all = TRUE) %>%
      arrange(desc(stat))
    ranks <- setNames(rank_df$stat, rank_df$gene)
    write.csv(rank_df %>% select(gene, stat, log2FoldChange, pvalue, padj),
              file.path(cmp_dir, "GSEA_rank_DESeq2_Wald_stat.csv"), row.names = FALSE)
    cat("Genes entering GSEA:", length(ranks), "\n")

    go_res <- run_gsea(ranks, go_pathways, cmp, "GO_BP", gsea_min_size[[pop_name]])
    re_res <- run_gsea(ranks, reactome_pathways, cmp, "Reactome", gsea_min_size[[pop_name]])

    write.csv(go_res, file.path(cmp_dir, "GO_BP_GSEA_lowExprFiltered_DESeq2.csv"), row.names = FALSE)
    write.csv(re_res, file.path(cmp_dir, "REACTOME_GSEA_lowExprFiltered_DESeq2.csv"), row.names = FALSE)
    write.csv(filter(go_res, FDR < 0.25, abs(NES) > 1.2),
              file.path(cmp_dir, "GO_BP_FDR025_NES12.csv"), row.names = FALSE)
    write.csv(filter(re_res, FDR < 0.25, abs(NES) > 1.2),
              file.path(cmp_dir, "REACTOME_FDR025_NES12.csv"), row.names = FALSE)

    gsea_results[[paste0(cmp, "_GO")]]       <- go_res
    gsea_results[[paste0(cmp, "_REACTOME")]] <- re_res
  }

  # ----------------------------------------------------------
  # 9.6 GSEA summary
  # ----------------------------------------------------------
  summary_tbl <- bind_rows(lapply(names(gsea_results), function(nm) {
    x   <- gsea_results[[nm]]
    cmp <- sub("_(GO|REACTOME)$", "", nm)
    tibble(
      population = pop_name,
      analysis   = nm,
      DEG_up     = deg_counts[[cmp]][["up"]],
      DEG_down   = deg_counts[[cmp]][["down"]],
      n_pathways = nrow(x),
      FDR005       = sum(x$FDR < 0.05, na.rm = TRUE),
      FDR025       = sum(x$FDR < 0.25, na.rm = TRUE),
      FDR025_NES12 = sum(x$FDR < 0.25 & abs(x$NES) > 1.2, na.rm = TRUE),
      positive_FDR025_NES12 = sum(x$FDR < 0.25 & x$NES >  1.2, na.rm = TRUE),
      negative_FDR025_NES12 = sum(x$FDR < 0.25 & x$NES < -1.2, na.rm = TRUE)
    )
  }))
  write.csv(summary_tbl, file.path(out_dir, "GSEA_summary.csv"), row.names = FALSE)
  print(summary_tbl, n = Inf)

  summary_tbl
}


# ============================================================
# 10. RUN ALL POPULATIONS
# ============================================================

all_summary <- bind_rows(lapply(names(populations), function(p) {
  run_population(p, populations[[p]])
}))


# ============================================================
# 11. MASTER SUMMARY
# ============================================================

cat("\n========================================\n")
cat("MASTER SUMMARY (all populations)\n")
cat("========================================\n")
print(all_summary, n = Inf)

write.csv(all_summary, file.path(OUT_ROOT, "GSEA_master_summary.csv"), row.names = FALSE)


# ============================================================
# 12. SESSION INFO
# ============================================================

writeLines(capture.output(sessionInfo()), file.path(OUT_ROOT, "sessionInfo.txt"))

cat("\n\nDONE.\nAll results saved to:\n", OUT_ROOT, "\n")
