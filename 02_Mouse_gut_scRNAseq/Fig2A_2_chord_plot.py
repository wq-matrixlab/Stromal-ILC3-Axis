import csv
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import ktplotspy as kpy

def setup_plot_environment():
   
    sc.settings.verbosity = 3
    sc.settings.set_figure_params(
        dpi=80,
        dpi_save=300,
        fontsize=14,
        format='pdf',
        facecolor='none'
    )

   
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']

def load_and_preprocess_data(means_path, pvals_path):
  
    means = pd.read_csv(means_path)
    pvals = pd.read_csv(pvals_path)

    numeric_means_cols = means.select_dtypes(include=[np.number, 'float', 'int']).columns
    for col in numeric_means_cols:
        means[col] = pd.to_numeric(means[col], errors='coerce')
        means[col] = means[col].where(means[col] >= 0.05, 0)

  
    numeric_pvals_cols = pvals.select_dtypes(include=[np.number, 'float', 'int']).columns
    for col in numeric_pvals_cols:
        pvals[col] = pd.to_numeric(pvals[col], errors='coerce')
        pvals[col] = pvals[col].where(pvals[col] <= 0.0001, 1)
        
    return means, pvals

def load_edge_colors(csv_path):
   
    color_map = {}
    with open(csv_path, 'r', encoding='utf-8-sig') as file:
        reader = csv.reader(file, delimiter=',')  
        for row in reader:
            if len(row) == 2:  
                key, value = row[0].strip(), row[1].strip()
                color_map[key] = value
    return color_map

def main():
  
    MEANS_FILE = 'means.csv'
    PVALS_FILE = 'pvalues.csv'
    DECON_FILE = 'statistical_analysis_deconvoluted.txt'
    ADATA_FILE = 'chord.h5ad'
    EDGE_COLOR_FILE = 'color_map.csv'
    OUTPUT_FILE = 'chord.pdf'

   
    FACE_COLORS = {
        'Monocyte': 'cornflowerblue',
        'CD4+T': 'pink',
        'CD8+T': 'palegreen',
        'ILC2': 'mediumslateblue',
        'ILC3': 'orchid',
        'Fibroblast(C3+)': 'lightcoral',
        'B cell': 'orange',
        'Neutrophil': 'coral',
        'Plasma cell': 'burlywood',
        'Macrophage': 'paleturquoise',
        'NK': 'turquoise',
        'DC': 'yellowgreen',
    }


    setup_plot_environment()
    
    print("Loading datasets...")
    means, pvals = load_and_preprocess_data(MEANS_FILE, PVALS_FILE)
    decon = pd.read_table(DECON_FILE)
    adata = sc.read_h5ad(ADATA_FILE)
    edge_color_map = load_edge_colors(EDGE_COLOR_FILE)

  
    print("\nProcessed Means Head:")
    print(means.head())
    print("\nProcessed P-values Head:")
    print(pvals.head())
    print("\nEdge Color Map Configuration loaded successfully.")

  
    print("\nGenerating chord diagram...")
    p = kpy.plot_cpdb_chord(
        adata=adata,
        cell_type1="Fibroblast(C3+)",  
        cell_type2=".",               
        means=means,
        pvals=pvals,
        deconvoluted=decon,
        celltype_key="cell type",
        face_col_dict=FACE_COLORS,
        edge_col_dict=edge_color_map,
        remove_self=True,             
        gap=0,
        scale_lw=3.5,                 
        size=10, 
        interspace=0,
        raxis_range=(940, 1000),       
        label_visible=True, 
        figsize=(5, 5), 
        legend_params={
            'bbox_to_anchor': (1, 1), 
            'frameon': False, 
            'loc': 'center left',
            'ncol': 1
        }, 
        layer=None, 
        alpha=0.5,
    ) 

  
    p.save(OUTPUT_FILE)
    print(f"Chord diagram successfully saved to {OUTPUT_FILE}")

if __name__ == '__main__':
    main()