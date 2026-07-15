import numpy as np  
import pandas as pd
import scanpy as sc
import anndata as ad
import scvelo as scv
import seaborn as sb
import matplotlib.pyplot as plt 

adata_annotation=sc.read_h5ad("Fibroblast_annotation.h5ad")
adata_velo=sc.read_h5ad("Fibroblast_velocity_counts.h5ad")
adata = scv.utils.merge(adata_annotation, adata_velo)

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=1000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.recover_dynamics(adata, n_jobs=16)

scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap', color='cell type', legend_loc='right margin')

ax = scv.pl.velocity_embedding_stream(
    adata,
    basis='umap',
    color='cell type',
    legend_loc='right margin',
    show=False
)

ax.set_axis_off()

plt.savefig("fibro_velo.png", bbox_inches="tight")
plt.show()

adata.write(
    "Fibroblast_scVelo_dynamical.h5ad"
)