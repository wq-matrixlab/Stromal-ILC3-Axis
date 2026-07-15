import scanpy as sc


INPUT_FILE = "UC_human_integrated_raw_counts.h5ad"
OUTPUT_FILE = "UC_human_QC_filtered.h5ad"



def calculate_qc_metrics(adata):




    adata.var["MT"] = (
        adata.var_names
        .str.upper()
        .str.startswith("MT-")
    )


    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["MT"],
        percent_top=None,
        log1p=False,
        inplace=True
    )


    return adata



def visualize_qc_metrics(adata):




    sc.pl.violin(
        adata,
        keys=[
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_MT"
        ],
        jitter=0.4,
        multi_panel=True
    )


    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_MT"
    )


    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts"
    )



def filter_cells_by_qc(
    adata,
    min_genes=500,
    min_cells=3,
    max_genes=6000,
    max_counts=60000,
    max_mt=20
):



    initial_cells = adata.n_obs
    initial_genes = adata.n_vars



    sc.pp.filter_cells(
        adata,
        min_genes=min_genes
    )



    sc.pp.filter_genes(
        adata,
        min_cells=min_cells
    )


    adata = adata[
        adata.obs["n_genes_by_counts"] < max_genes
    ].copy()



    adata = adata[
        adata.obs["total_counts"] < max_counts
    ].copy()



    adata = adata[
        adata.obs["pct_counts_MT"] < max_mt
    ].copy()


    print(
        "QC filtering summary:"
    )

    print(
        f"Cells: {initial_cells} → {adata.n_obs}"
    )

    print(
        f"Genes: {initial_genes} → {adata.n_vars}"
    )


    return adata



def main():

    adata = sc.read_h5ad(
        INPUT_FILE
    )


    print(
        f"Initial dataset: {adata.shape}"
    )


    adata = calculate_qc_metrics(
        adata
    )


    visualize_qc_metrics(
        adata
    )


    adata = filter_cells_by_qc(
        adata
    )


    adata.write_h5ad(
        OUTPUT_FILE,
        compression="gzip"
    )


if __name__ == "__main__":

    main()