# MASH atrial snRNA-seq

Code for the single-nucleus RNA-seq analysis accompanying the manuscript.
We sequenced atria from
Chow, MASH, MASH + anti-OPN, MASH *Spp1*<sup>f/f</sup> and MASH
*Spp1*<sup>f/f</sup>;*Alb*<sup>Cre</sup> mice (13 libraries in total).

Raw and processed data are in GEO under GSE331217.

## What is here

`scripts/` has the analysis code, numbered in the order we ran it. metadata/ has the sample sheet (`sample_metadata.tsv`), final cell annotations
(`cell_annotations.csv.gz`), MP/DC subtype labels (`mpdc_subtypes.csv.gz`) and
fibroblast subtype labels (`fb_subtypes.csv.gz`). `logs/` keeps the CellBender parameters,
the per-sample QC summary and the number of nuclei left after each filtering step.

## Running it

Each sample needs its own folder with the Cell Ranger output, and this repository is
expected inside the same project folder (`<project>/MASH-snRNAseq`). Set `DATA_DIR` at
the top of the R and Python scripts to the project folder. The CellBender script reads the
same path from `MASH_DATA_DIR`:

```bash
export MASH_DATA_DIR=/path/to/project
```

Intermediate `.rds`/`.h5ad` files go to the project folder and figures/tables go to
`<project>/results/`. None of these are tracked here.

## Scripts

Preprocessing and QC

- `00b_cellbender_commands.sh` – CellBender 0.3.0, run per sample with the settings we used
- `00c_cellbender_h5_to_h5ad.py` – converts the CellBender output to h5ad
- `01_QC_doublets.R` – per-sample filtering and scDblFinder
- `02_merge.R` – merges the 13 samples
- `03_export_for_scvi.R`, `04_scvi_integration.py` – scVI integration (batch = sample)
- `05_QC_annotation.R` – clustering, annotation and removal of contaminating or
  low-quality nuclei; writes `seurat_object_QC4.rds` and `cell_annotations.csv.gz`

Analysis

- `06_overview_figures.R` – Harmony embedding of the final object (`seurat_final.rds`),
  UMAP, cell-type proportions, marker dot plot
- `07_pseudobulk_PCA.R` – sample-level PCA
- `08_MiloR.R` – differential abundance (miloR, `~ batch + group`)
- `09a_MPDC_subclusters.R` – MP1–MP4, Mono and DC; also exports the table used for velocity
- `09b_RNA_velocity_MPDC.py` – scVelo on monocytes and macrophages
- `10_fibroblast_subclusters.R` – FB1–FB5, activation score, module scores,
  pseudobulk DEGs
- `11_CellChat.R` – signalling between MP1, FB1, atrial cardiomyocytes and endothelial cells
- `12_pseudobulk_DESeq2_GSEA_all_populations.R` – pseudobulk DESeq2 and GSEA for whole
  atrium, FB, MP/DC, EC and ACM
- `13_final_figures.R` – builds the manuscript figures from the outputs above

11 needs the objects from 09a and 10. 13 needs 06, 09a and 12.

## Notes on the methods

Pseudobulk counts are summed per sample. Genes were filtered with `edgeR::filterByExpr`
(min.count = 10, min.total.count = 15) and tested in DESeq2 with all 13 samples and
`~ batch + group`. Anti-OPN mice were partly processed in their own batch, so the
anti-OPN vs MASH contrast should be read with that in mind.

For GSEA, genes were ranked by the DESeq2 Wald statistic and run through
fgseaMultilevel against MSigDB GO:BP and Reactome (human sets mapped to mouse
orthologs). maxSize was 500; minSize was 10 for whole atrium and MP/DC and 15 for FB, EC
and ACM, matching the runs behind the figures. Mitochondrial, ribosomal, hemoglobin and
sex-linked genes were left out of the ranking but not out of DESeq2.

For RNA velocity, loom files were made with velocyto from the Cell Ranger BAMs; looms
from multiplexed runs were split using each sample's CellBender barcodes. That step is
not included. Spliced/unspliced counts were not batch-corrected, moments were computed on
a PCA kNN graph, and the vectors are drawn on the Harmony UMAP.

The clusters dropped in 05 were picked by looking at marker genes. Cluster numbers can
shift with other package versions, so `cell_annotations.csv.gz` is the reference for
which nuclei and labels were used in the paper.

## Software

R 4.x with Seurat 5, harmony, scDblFinder, zellkonverter, limma, edgeR, DESeq2, fgsea,
msigdbr, CellChat 2, miloR, Nebulosa, ComplexHeatmap, circlize and ggnewscale.
Python 3 with CellBender 0.3.0, scanpy, scvi-tools and scVelo.
Exact Python package versions are provided in `logs/python_versions.txt`.
R package and session information is provided in `logs/sessionInfo_R.txt`.

## Contact

Code questions: Jifei Ding (Jifei.Ding@bcm.edu), or open an issue on this repository.

Corresponding author: Na Li (nal@bcm.edu)
