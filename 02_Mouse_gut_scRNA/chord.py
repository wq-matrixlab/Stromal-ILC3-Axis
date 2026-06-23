import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import ktplotspy as kpy

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
means = pd.read_csv('means.csv')
pvals = pd.read_csv('pvalues.csv')
decon = pd.read_table("statistical_analysis_deconvoluted.txt")
adata = sc.read_h5ad('chord.h5ad')

for col in means.columns:
    if pd.api.types.is_numeric_dtype(means[col]):
        means[col] = pd.to_numeric(means[col], errors='coerce')
        means[col] = means[col].where(means[col] >= 0.05, 0)

print(means.head())

for col in pvals.columns:
    if pd.api.types.is_numeric_dtype(pvals[col]):
        pvals[col] = pd.to_numeric(pvals[col], errors='coerce')
        pvals[col] = pvals[col].where(pvals[col] <= 0.0001, 1)

print(pvals.head())
import csv

color_map = {}


with open('AD5.csv', 'r', encoding='utf-8-sig') as file:
    reader = csv.reader(file, delimiter=',')  
    for row in reader:
        if len(row) == 2:  
            key = row[0]  
            value = row[1]  
            color_map[key] = value

print(color_map)
p = kpy.plot_cpdb_chord(
    adata=adata,
    cell_type1="Fibroblast(C3+)",
    cell_type2=".",
    means=means,
    pvals=pvals,
    deconvoluted=decon,
    celltype_key="cell type",
    edge_cmap=plt.cm.coolwarm,
    face_col_dict = {
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
       } ,
    remove_self=True,
    gap=0,
    scale_lw=3.5,
    size=10, 
    interspace=0,
    raxis_range=(940, 1000),
  
    label_visible=True, 
    figsize=(5, 5), 
    legend_params={'bbox_to_anchor': (1, 1), 
                   'frameon': False, 'loc': 'center left','ncol':1}, 
    layer=None, 
    edge_col_dict=color_map,
  
    alpha=0.5,
    ) 
p.save('chord.pdf') 