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
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.set_figure_params(dpi=50, fontsize=12, dpi_save=300)

fig_outdir= '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_fig'
directory = '/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas'
os.chdir(directory)
datadir = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scVI/'

adata= sc.read_h5ad(os.path.join(datadir,'scVI_adata_anno.h5ad'))

adata_singlet = adata[adata.obs['doublet'] == 'singlet'].copy()