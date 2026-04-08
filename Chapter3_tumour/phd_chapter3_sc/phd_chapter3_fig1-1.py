import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns 
import matplotlib.pyplot as plt
import scipy
import anndata as ad
from pathlib import Path
import glob
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec

sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.set_figure_params(dpi=300, fontsize=12, dpi_save=300)

from matplotlib.image import imread

# Path
fig_outdir= '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_fig'
directory = '/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas'
os.chdir(directory)
datadir = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scVI/'

adata_or= sc.read_h5ad(os.path.join(datadir,'scVI_adata_anno.h5ad'))

adata_singlet = adata_or[adata_or.obs['doublet'] == 'singlet'].copy()
pd.crosstab(adata_singlet.obs['source'],adata_singlet.obs['custom_annotation'])

# Color dictionary
color_dict = {
    'Tumor': '#ff7f0eff',
    'OPCs': '#ff0000ff',
    'Neurons': '#17becfff',
    'Immune cells': '#aa40fcff',
    'Oligodendrocytes': '#279e68ff'
}

categories = adata_singlet.obs['custom_annotation'].cat.categories

# Assign colors to .uns
adata_singlet.uns['custom_annotation_colors'] = [color_dict[cat] for cat in categories]

# Plot per cell type
t= sc.pl.umap(
    adata_singlet,
    color='custom_annotation',
    frameon=False,
    title='Cell type'
)


ct = pd.crosstab(adata_singlet.obs["source"],adata_singlet.obs["custom_annotation"])
ct_props = ct.div(ct.sum(axis=1), axis=0)

# Normalize counts to proportions
ct_props = ct.div(ct.sum(axis=1), axis=0)

# Plot
fig, ax = plt.subplots(figsize=(5, 3), facecolor='white')

ct_props.plot(
    kind='bar',
    stacked=True,
    color=[color_dict[col] for col in ct.columns],
    ax=ax,
    width=0.95,          # Thin space between bars
    edgecolor='white'    # White separator between segments
)

# Style
ax.set_ylim(0, 1)
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_yticklabels(['0', '25', '50', '75', '100'])
ax.set_xlabel('')
for i, label in enumerate(ax.get_xticklabels()):
    if i % 2 == 0:
        label.set_fontweight('normal')  # odd index (0-based) labels normal
    else:
        label.set_fontweight('bold') 
ax.set_title('Cell Type Proportions (%)')
ax.grid(False)
ax.legend_.remove()  # Remove legend

# Remove border frame
for spine in ax.spines.values():
    spine.set_visible(False)

fig.savefig(os.path.join(fig_outdir,"cell_type_proportions.png"), dpi=300, bbox_inches='tight', transparent=True)

plt.tight_layout()
plt.show()

# Dot plot
sc.pl.dotplot(
    adata_singlet,
    ['MOG', 'MBP', 'PLP1', 'SNAP25', 'RBFOX1', 'RBFOX3', 'DOCK8', 'PTPRC', 'MERTK', 'CD74', 'CD163'],
    groupby="custom_annotation",
    categories_order=['Immune cells', 'Neurons', 'Oligodendrocytes', 'Transitional', 'Tumor'])


# SM fig1: Plot 4 panels
ax = sc.pl.umap(
    adata_or,
    color=['source','Pair', 'Tumor_Normal_annotation','Status'],
    frameon=False,
    ncols=2,
    wspace=0.5,
    show=False
)

# Access figure from axes
fig = ax[0].figure

# Custom panel titles
titles = ['Sample ID','Patient pair','Malignant status', 'Recurent state']

# Assign titles to each subplot
for a, title in zip(ax, titles):
    a.set_title(title, fontsize=14, weight='bold')

fig.savefig(os.path.join(fig_outdir, 'gbm_panels_labeled.png'), dpi=300, bbox_inches='tight')
