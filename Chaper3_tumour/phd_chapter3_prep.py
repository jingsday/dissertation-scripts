import os
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns 
import matplotlib.pyplot as plt
import scipy
import csv
import gzip
import anndata as ad
from pathlib import Path
import glob
import scvi
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import ListedColormap
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

directory = '/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas'
os.chdir(directory)
outdir = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_publication/'

names_list=['GSM5319518_SF2777','GSM5319548_SF2979','GSM5319519_SF2990',
                'GSM5319549_SF3073','GSM5319520_SF3076','GSM5319550_SF3243',
                'GSM5319521_SF3391','GSM5319551_SF3448','GSM5319511_SF11916',
                'GSM5319543_SF12382','GSM5319506_SF11082','GSM5319562_SF11488',
                'GSM5319530_SF9358','GSM5319568_SF9962','GSM5319559_SF9798','GSM5319532_SF9494']

adata_list = []

# Loop over each sample and read in the AnnData object
for name in names_list:

    mtx =f"{name}_matrix.mtx.gz"
    adata = sc.read_mtx(mtx)
    cells=pd.read_csv(f'{name}_barcodes.tsv.gz',header=None)
    features=pd.read_csv(f'{name}_features.tsv.gz',header=None,sep='\t')
    adata= adata.T
    #check the columns first to make sure they are the ones you need 
    adata.obs['CellID']= cells[0].tolist()
    adata.obs.index = adata.obs['CellID']
    adata.var['Gene']= features[0].tolist()
    adata.var.index= adata.var['Gene']
    adata.var_names_make_unique() 
    adata.var['mt'] =adata.var_names.str.startswith('MT-')


    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata= adata[adata.obs.n_genes_by_counts <6000, :]
    adata= adata[adata.obs.pct_counts_mt< 5, :].copy()


    adata.obs['source'] = name[11:]
    adata.layers["counts"] = adata.X.copy()    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep full dimension safe
    adata_list.append(adata)

batch_names = [adata.obs['source'].iloc[0] for adata in adata_list]
adata = adata_list[0].concatenate(adata_list[1:], batch_key='source', batch_categories=batch_names)

adata.raw = adata  # keep full dimension safe
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="source",
    subset=True,
)

# scVI 
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="source")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train()

SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()

adata.obs = adata.obs.rename(columns={'CellID': 'Barcodes'})

adata.write_h5ad(os.path.join(outdir,'scVI_adata.h5ad'),compression='gzip')

sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata,resolution=0.5)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color=["source", "leiden"],
    frameon=False,
    ncols=1,
)

model_outdir = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUTPUT_scVI/'
model_dir = os.path.join(model_outdir, "scvi_model")
model.save(model_dir, overwrite=True)
#scvi.model.SCVI.load(model_dir, adata=adata)

# Doublet identification
doublet_df = {}  

for i in set(adata.obs['source']):
    d = scvi.external.SOLO.from_scvi_model(model, restrict_to_batch=i)
    d.train()
    
    df_bulk = d.predict()
    
    df_bulk['solo_prediction'] = d.predict(soft=False)
    df_bulk.index = df_bulk.index.map(lambda x: x)
    
    doublet_df.update(dict(zip(df_bulk.index, df_bulk['solo_prediction'])))

df_output = pd.DataFrame(list(doublet_df.items()), columns=['Index', 'Solo_Prediction'])
df_output.to_csv(os.path.join(outdir,'doublet_predictions.csv'), index=False)

for i in adata.obs.index:
    adata.obs.loc[i,'doublet'] = doublet_df.loc[i,'solo']

# Metadata
patients_info = pd.read_excel('/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_sm.xlsx',sheet_name='Table 1',index_col=0)
patients_info
patients_info.reset_index(inplace=True)
adata.obs.reset_index(inplace= True)
adata.obs= adata.obs.merge(patients_info[['ID','Stage', 'Pair#', 'Sex', 'Age']], right_on='ID',left_on='source',how='left')

adata.obs.set_index('CellID',inplace=True)

adata.obs["Pair#"] = adata.obs["Pair#"].astype(str)
adata.write_h5ad(os.path.join(outdir, 'scVI_adata.h5ad'), compression='gzip')
