import os
import numpy as np
import scanpy as sc
import anndata as ad

from scipy import sparse


DATA_DIR = "sc_UC_human"

DATASETS = [
    "Clean_JK_Cell_2018_raw_count.h5ad",
    "Clean_C_S_Smillie_et_al_Cell_2020.h5ad",
    "Clean_M_Friedrich_Nat_Med_2021.h5ad",
    "Clean_Rasa_Elmentaite_et_al_Nature_2021.h5ad"
]

OUTPUT = "UC_human_integrated_raw_counts.h5ad"



def load_datasets():

    adatas=[]

    for file in DATASETS:

        path=os.path.join(DATA_DIR,file)

        if not os.path.exists(path):
            raise FileNotFoundError(path)

        print(f"Loading {file}")

        adatas.append(
            sc.read_h5ad(path)
        )

    return ad.concat(
        adatas,
        join="outer",
        merge="same",
        uns_merge="same"
    )



def clean_adata(adata):

    print(
        f"Before filtering: {adata.shape}"
    )


    if sparse.issparse(adata.X):

        adata.X.data=np.nan_to_num(
            adata.X.data,
            nan=0
        )

    else:

        adata.X=np.nan_to_num(
            adata.X,
            nan=0
        )


    keep_obs=[
        c for c in 
        ['class','Sample','RawData']
        if c in adata.obs.columns
    ]

    adata.obs=adata.obs[keep_obs]


    adata=adata[
        adata.obs['class']
        .isin(
            ['Healthy','UC_inflamed']
        )
    ]


    print(
        f"After filtering: {adata.shape}"
    )

    return adata



if __name__=="__main__":

    adata=load_datasets()

    adata=clean_adata(adata)

    adata.write_h5ad(
        OUTPUT,
        compression="gzip"
    )

    print(
        f"Saved: {OUTPUT}"
    )