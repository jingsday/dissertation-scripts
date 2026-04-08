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

adata= sc.read_h5ad('/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_publication/scVI_adata_anno.h5ad')
adata_singlet = adata[adata.obs['solo']=='singlet',].copy()

adata_singlet.obs['DEG_anno'] = adata_singlet.obs['Status'].astype(str)+'_'+ adata_singlet.obs['custom_annotation'].astype(str)

sc.tl.rank_genes_groups(
    adata_singlet,
    groupby="DEG_anno",        
    method="wilcoxon",       
    use_raw=True             
)

sc.tl.rank_genes_groups(adata_singlet, "DEG_anno", groups=["Recurrent_Oligodendrocytes"], reference="Primary_Oligodendrocytes", method="wilcoxon")
sc.pl.rank_genes_groups(adata_singlet, groups=["Recurrent_Oligodendrocytes"], n_genes=20)

result = adata_singlet.uns["rank_genes_groups"]
groups = result["names"].dtype.names
df= pd.DataFrame(
    {
        f"{group}_{key[:1]}": result[key][group]
        for group in groups
        for key in ["names",'logfoldchanges', "pvals"]
    }
)   

df.columns =['names', 'logfoldchanges','pvals']



# rhoa genes scores
rhoa_genes = [
    "RHOA","RHOB","RHOC","ROCK1","ROCK2",
    "MYL9","MYH9","PPP1R12A","DIAPH1","PFN1","ACTN1","TAGLN","CNN1"
]

sc.tl.score_genes(adata_singlet, gene_list=rhoa_genes, score_name='RhoA_score')

import pandas as pd
from scipy.stats import ttest_rel

# Subset neurons
neurons = adata_singlet[adata_singlet.obs['custom_annotation'] == 'Neurons'].obs.copy()

# Clean labels
neurons['Pair'] = neurons['Pair'].str.strip()
neurons['Stage'] = neurons['Stage'].str.lower().str.strip()

# Compute per-patient mean RhoA scores
mean_scores = neurons.groupby(['Pair', 'Stage'])['RhoA_score'].mean().unstack(fill_value=pd.NA)

print(mean_scores)  # check what you have

# Keep only patients that have both primary and recurrent
paired_scores = mean_scores.dropna(subset=['primary', 'recurrent'])

# Run paired t-test
t_stat, p_val = ttest_rel(paired_scores['primary'], paired_scores['recurrent'])
print("t =", t_stat, "p =", p_val)
