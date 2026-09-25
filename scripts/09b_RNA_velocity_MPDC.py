"""
RNA velocity (scVelo, stochastic model) of monocytes and macrophage sub-types.

Input : <VEL_DIR>/looms/<Sample>.loom
        Pre-generated sample-level velocyto loom files
        <VEL_DIR>/metadata_MP_velocity.csv
"""
import csv
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

DATA_DIR = "/path/to/MASHproject"
REPO_DIR = os.path.join(DATA_DIR, "MASH-snRNAseq")
VEL_DIR = os.path.join(DATA_DIR, "velocity")
FIG_DIR = os.path.join(VEL_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

scv.settings.verbosity = 2
scv.settings.set_figure_params("scvelo", dpi=300)

COLORS = {"Mono": "#7570B3", "MP1": "#E7298A", "MP2": "#66A61E",
          "MP3": "#E6AB02", "MP4": "#1B9E77"}
PANELS = {"CHOW": "Chow", "MASH": "MASH",
          "MASH_flox": "Spp1 f/f", "MASH_SPP1KO": "Spp1 f/f;Alb-Cre"}


def save(name):
    plt.savefig(os.path.join(FIG_DIR, name + ".png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(FIG_DIR, name + ".svg"), bbox_inches="tight")
    plt.close()


def run_velocity(adata):
    """Normalisation, PCA-based kNN graph, moments, stochastic velocity."""
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)        # kNN graph on X_pca
    scv.tl.velocity(adata, mode="stochastic")
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_embedding(adata, basis="umap")
    return adata


def plot_group_arrows(adata, prefix):
    for grp, title in PANELS.items():
        sub = adata[adata.obs["group"] == grp].copy()   # keeps velocity_umap of the joint fit
        if sub.n_obs == 0:
            continue
        scv.pl.velocity_embedding(sub, basis="umap", color="celltype",
                                  arrow_length=3, arrow_size=2, title=title, show=False)
        save(f"{prefix}_{grp}")


# -----------------------------------------------------------------------------
# 1. Read per-sample looms; cell names = <orig.ident>_<BARCODE>
# -----------------------------------------------------------------------------
with open(os.path.join(REPO_DIR, "metadata", "sample_metadata.tsv")) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
samples = {r["Sample"]: r["orig.ident"] for r in rows}
batch_of = {r["orig.ident"]: r["Batch"] for r in rows}


def clean_barcode(cell_id, sample):
    bc = str(cell_id).split(":")[-1].rstrip("x")      # velocyto: <run>:<BARCODE>x
    if bc.startswith(sample + "_"):                     # split looms: <Sample>_<BARCODE>
        bc = bc[len(sample) + 1:]
    return bc.replace("-1", "")


adatas = []
for sample, oi in samples.items():
    a = sc.read_loom(os.path.join(VEL_DIR, "looms", f"{sample}.loom"), sparse=True)
    a.var_names_make_unique()
    a.obs_names = [f"{oi}_{clean_barcode(c, sample)}" for c in a.obs_names]
    adatas.append(a)
adata = ad.concat(adatas, merge="same")

# -----------------------------------------------------------------------------
# 2. Keep monocytes/macrophages annotated in Seurat; add labels and UMAP
# -----------------------------------------------------------------------------
meta = pd.read_csv(os.path.join(VEL_DIR, "metadata_MP_velocity.csv"))
meta.index = meta["cell_id"].str.replace("-1$", "", regex=True)
common = meta.index.intersection(adata.obs_names)
print(f"Matched {len(common)} / {len(meta)} cells")
adata = adata[common].copy()
meta = meta.loc[common]

for col in ["orig.ident", "group", "celltype"]:
    adata.obs[col] = meta[col].values
adata.obs["batch"] = adata.obs["orig.ident"].map(batch_of).astype("category")
adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"], categories=list(COLORS))
adata.uns["celltype_colors"] = list(COLORS.values())
adata.obsm["X_umap"] = meta[["UMAP_1", "UMAP_2"]].to_numpy()   # Harmony-based MP/DC UMAP
raw = adata.copy()                                            # for the batch3 refit

# -----------------------------------------------------------------------------
# 3. Velocity on all MP/DC nuclei
# -----------------------------------------------------------------------------
adata = run_velocity(adata)
scv.tl.velocity_confidence(adata)
adata.write(os.path.join(VEL_DIR, "MPDC_velocity.h5ad"))

# -----------------------------------------------------------------------------
# 4. Diagnostics: is the PCA kNN graph driven by batch?
# -----------------------------------------------------------------------------
scv.pl.scatter(adata, basis="pca", color=["batch", "celltype", "group"], ncols=3, show=False)
save("diagnostic_PCA_batch_celltype_group")

knn = adata.obsp["distances"].tocsr()
b = adata.obs["batch"].to_numpy()
same = np.mean([np.mean(b[knn[i].indices] == b[i]) for i in range(adata.n_obs)])
expected = (adata.obs["batch"].value_counts(normalize=True) ** 2).sum()
print(f"Same-batch neighbours: {same:.2f} (random expectation {expected:.2f})")

# -----------------------------------------------------------------------------
# 5. Figures
# -----------------------------------------------------------------------------
scv.pl.velocity_embedding_stream(adata, basis="umap", color="celltype",
                                 title="RNA velocity", show=False)
save("velocity_stream_all")

scv.pl.scatter(adata, basis="umap", color=["velocity_length", "velocity_confidence"],
               color_map="coolwarm", show=False)
save("velocity_confidence")

plot_group_arrows(adata, "velocity_arrow")

# -----------------------------------------------------------------------------
# 6. Robustness: refit within batch3 (contains all five groups)
# -----------------------------------------------------------------------------
b3 = run_velocity(raw[raw.obs["batch"] == "batch3"].copy())
scv.pl.velocity_embedding_stream(b3, basis="umap", color="celltype",
                                 title="RNA velocity (batch3 only)", show=False)
save("robustness_batch3_stream")
plot_group_arrows(b3, "robustness_batch3_arrow")
