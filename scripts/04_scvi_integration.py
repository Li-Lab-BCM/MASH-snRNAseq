"""
04_scvi_integration.py
scVI integration across samples (batch = orig.ident).

Input : <DATA_DIR>/seurat_for_scvi_QC1.h5ad   (raw counts)
Output: <DATA_DIR>/seurat_scvi_5group.h5ad    (obsm["X_scVI"], 30 dims)
        <DATA_DIR>/scvi_model/
"""
import os

import anndata
import scanpy as sc
import scvi

DATA_DIR = "/path/to/MASHproject"
scvi.settings.seed = 0

adata = anndata.read_h5ad(os.path.join(DATA_DIR, "seurat_for_scvi_QC1.h5ad"))
if "counts" in adata.layers:
    adata.X = adata.layers["counts"]
print("cells:", adata.n_obs, "genes:", adata.n_vars)

# Highly variable genes (on raw counts, per-sample aware)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="orig.ident",
    flavor="seurat_v3",
    subset=True,
)

scvi.model.SCVI.setup_anndata(adata, batch_key="orig.ident")
model = scvi.model.SCVI(
    adata,
    n_layers=2,
    n_latent=30,
    n_hidden=512,
    gene_likelihood="nb",
    dropout_rate=0.1,
)
model.train(
    max_epochs=400,
    batch_size=1024,
    early_stopping=True,
    early_stopping_patience=20,
    plan_kwargs={"lr": 1e-3},
)

adata.obsm["X_scVI"] = model.get_latent_representation()
adata.write(os.path.join(DATA_DIR, "seurat_scvi_5group.h5ad"))
model.save(os.path.join(DATA_DIR, "scvi_model"), overwrite=True)
print("Done:", adata.obsm["X_scVI"].shape)
