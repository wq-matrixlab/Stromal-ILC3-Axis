import scanpy as sc
import harmonypy as hm
import numpy as np


INPUT_FILE = "UC_human_QC_filtered.h5ad"
OUTPUT_FILE = "Fibroblast_PCA_harmony.h5ad"


sc.settings.verbosity = 3

sc.settings.set_figure_params(
    dpi=80,
    dpi_save=300,
    fontsize=8,
    figsize=(4,4),
    facecolor="white",
    transparent=True
)



adata = sc.read_h5ad(
    INPUT_FILE
)


print(
    f"Loaded data: {adata.shape}"
)



sc.pp.normalize_total(
    adata,
    target_sum=1e4
)


sc.pp.log1p(
    adata
)


sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_top_genes=2000,
    batch_key="Sample"
)


sc.pl.highly_variable_genes(
    adata
)

adata.raw = adata.copy()


adata = adata[
    :,
    adata.var.highly_variable
].copy()


sc.pp.regress_out(
    adata,
    [
        "total_counts",
        "pct_counts_MT"
    ]
)


sc.pp.scale(
    adata,
    max_value=10
)


sc.tl.pca(
    adata,
    svd_solver="arpack"
)

harmony_out = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    ["Sample"]
)


adata.obsm["X_pca_harmony"] = (
    harmony_out.Z_corr.T
)

sc.pp.neighbors(
    adata,
    use_rep="X_pca_harmony",
    n_neighbors=50,
    n_pcs=40,
    random_state=7
)


sc.tl.umap(
    adata,
    random_state=7
)


sc.tl.leiden(
    adata,
    resolution=0.2,
    random_state=0
)

sc.pl.umap(
    adata,
    color="leiden"
)


sc.pl.umap(
    adata,
    color="leiden",
    legend_loc="on data"
)

marker_genes = [
    "CDH5",
    "PECAM1",
    "SOX10",
    "PLP1",
    "PDGFRA",
    "DCN",
    "LUM",
    "COL1A1",
    "COL3A1",
    "ACTA2",
    "TAGLN",
    "MYH11",
    "PDGFRB",
    "RGS5",
    "PLN",
    "EPCAM",
    "MUC2",
    "PTPRC"
]


sc.pl.dotplot(
    adata,
    marker_genes,
    groupby="leiden",
    dendrogram=True,
    swap_axes=True
)

adata = adata.raw.to_adata()


adata = adata[
    np.isin(
        adata.obs["leiden"],
        ["3"]
    )
].copy()


sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    n_top_genes=2000,
    batch_key="Sample"
)


sc.pl.highly_variable_genes(
    adata
)

adata.raw = adata.copy()


adata = adata[
    :,
    adata.var.highly_variable
].copy()


sc.pp.regress_out(
    adata,
    [
        "total_counts",
        "pct_counts_MT"
    ]
)


sc.pp.scale(
    adata,
    max_value=10
)


sc.tl.pca(
    adata,
    svd_solver="arpack"
)

harmony_out = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    ["Sample"],max_iter_harmony=50
)


adata.obsm["X_pca_harmony"] = (
    harmony_out.Z_corr.T
)

sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=50,n_pcs=30)
sc.tl.umap(adata, min_dist=0.2,spread=1)
sc.pl.umap(adata)

adata.write_h5ad(
    OUTPUT_FILE,
    compression="gzip"
)