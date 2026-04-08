# Environment pySCENIC 3.10.17
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp 

import json 
import zlib
import base64

# collect SCENIC AUCell output
DATASET_ID = 'gbm2'
RESULTS_FOLDERNAME = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUT_publication_scenic_FINAL'
f_pyscenic_output = os.path.join(RESULTS_FOLDERNAME, '{}_final_scenic.loom'.format(DATASET_ID))

lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)


from pyscenic.rss import regulon_specificity_scores

anno= sc.read_h5ad(os.path.join('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scVI/','scVI_adata_anno.h5ad'))

anno.obs['further_type'] = anno.obs['Status'].astype('str') + '_' + anno.obs['custom_annotation'].astype('str')
anno.obs['further_type'].value_counts()

rss_cellType = regulon_specificity_scores( auc_mtx, anno.obs['further_type'] )
rss_cellType.T['Primary_Oligodendrocytes'].sort_values(ascending=False).head(20).index
rss_cellType.T['Recurrent_Oligodendrocytes'].sort_values(ascending=False).head(20).index

import seaborn as sns
import matplotlib.pyplot as plt
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize

auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
#auc_mtx_Z.sort_index(inplace=True)


# Figure 1e
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 12,
    'font.family': 'DejaVu Sans',
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'axes.labelweight': 'bold',
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'legend.title_fontsize': 12,
})

clusters_to_plot = ['Primary_Oligodendrocytes','Recurrent_Oligodendrocytes']

cell_mask = anno.obs['further_type'].isin(clusters_to_plot)
filtered_cells = anno.obs[cell_mask]

# Select genes of interest
selected_genes = ['SOX17(+)', 'FOXL2(+)', 'ISL1(+)', 'TBX18(+)', 'PGR(+)', 'FEV(+)',
       'EVX2(+)', 'GATA3(+)', 'GFI1(+)', 'HOXA11(+)', 'IKZF2(+)', 'MYBL1(+)',
       'TBX21(+)', 'FOXL1(+)', 'SOX10(+)', 'STAT2(+)', 'SPIB(+)', 'ERG(+)',
       'MXI1(+)', 'WT1(+)',  'TFEB(+)', 'FOXN2(+)',
       'ELF1(+)', 'BACH2(+)', 'NFIX(+)', 'ELF2(+)', 'FOXP1(+)', 'MYRF(+)',
       'NFE2L2(+)', 'FOXO1(+)', 'GTF2I(+)', 'BPTF(+)', 'SOX2(+)', 'FLI1(+)',
       'IKZF1(+)', 'POU2F1(+)'] # duplicate 'SOX10(+)','MXI1(+)',IKZF2(+)','STAT2(+)', 


# Subset and sort AUC matrix
filtered_auc = auc_mtx_Z.loc[filtered_cells.index, selected_genes]
cluster_labels = filtered_cells['further_type']
sorted_index = cluster_labels.sort_values().index
sorted_auc = filtered_auc.loc[sorted_index]
sorted_cluster_labels = cluster_labels.loc[sorted_index]

# Define colors for clusters
unique_clusters = sorted_cluster_labels.unique()
palette = sns.color_palette("tab10", len(unique_clusters))
lut = dict(zip(unique_clusters, palette))
row_colors = sorted_cluster_labels.map(lut)

g = sns.clustermap(
    sorted_auc,
    row_cluster=True,
    col_cluster=True,
    yticklabels=False,
    xticklabels=True,
    row_colors=row_colors,
    linecolor='gray',
    cmap="vlag",
    vmin=-2,
    vmax=3,
    figsize=(24, 20),
    dendrogram_ratio=(.1, .2),
)
g.ax_heatmap.tick_params(axis='x', labelsize=20)  # or 16, 18, as you prefer
if g.cax is not None:
    g.cax.tick_params(labelsize=18)  # or any size you want
#g.cax.set_ylabel("Regulon activity z-score", fontsize=16, fontweight='bold')


g.fig.tight_layout()
plt.show()