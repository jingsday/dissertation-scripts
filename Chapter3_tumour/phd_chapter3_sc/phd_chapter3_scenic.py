import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

# 20 threads
sc.settings.njobs = 20

directory = '/home/jing/Phd_project/project_GBM/gbm_DATA/gbm_DATA_GSE174554/gbm_DATA_scRNA_atlas'
os.chdir(directory)

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
    adata_list.append(adata)

batch_names = [adata.obs['source'].iloc[0] for adata in adata_list]
adata = adata_list[0].concatenate(adata_list[1:], batch_key='source', batch_categories=batch_names) 
adata.obs.index.name = "ref"



wdir = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUT_publication_scenic_FINAL'
os.chdir( wdir )
# Target documents
f_loom_path_scenic = "gbm2_filtered_scenic.loom"
f_anndata_path = "anndata2.h5ad"
f_pyscenic_output = "pyscenic2_output.loom"
f_final_loom = 'gbm2_scenic_integrated-output.loom'

adata.write( f_anndata_path )

adata_singlet = adata[adata.obs['doublet']=='singlet',].copy()
adata_singlet

# loom file
row_attrs = {
    "Gene": np.array(adata_singlet.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata_singlet.obs_names) ,
    "nGene": np.array( np.sum(adata_singlet.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_singlet.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata_singlet.X.transpose(), row_attrs, col_attrs)


# SCENIC pipeline

# Resource
AUXILLIARIES_FOLDERNAME = "/home/jing/pySCENIC/resources/"
HUMAN_TFS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'allTFs_hg38.txt')

RANKING_DBS_FNAMES = list(map(lambda fn: os.path.join(AUXILLIARIES_FOLDERNAME, fn),
                       ['hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
                       'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather']))

MOTIF_ANNOTATIONS_FNAME = os.path.join(AUXILLIARIES_FOLDERNAME, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

#Outputs
DATASET_ID = 'gbm2'
RESULTS_FOLDERNAME = '/home/jing/Phd_project/project_GBM/gbm_OUTPUT/gbm_OUT_publication_scenic_FINAL'
os.chdir(RESULTS_FOLDERNAME)
F_LOOM_SCE =  os.path.join(RESULTS_FOLDERNAME, '{}_filtered_scenic.loom'.format(DATASET_ID))
f_pyscenic_output = os.path.join(RESULTS_FOLDERNAME, '{}_final_scenic.loom'.format(DATASET_ID))


ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(DATASET_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.motifs.csv'.format(DATASET_ID))

REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.regulons.dat'.format(DATASET_ID))
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc.csv'.format(DATASET_ID))
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.bin.csv'.format(DATASET_ID))
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.thresholds.csv'.format(DATASET_ID))

os.chdir(RESULTS_FOLDERNAME)

#I) GRN
!pyscenic grn {F_LOOM_SCE} {HUMAN_TFS_FNAME} -o {ADJACENCIES_FNAME} --num_workers 15

DBS_PARAM = ' '.join(RANKING_DBS_FNAMES)

#II) CTX
!pyscenic ctx {ADJACENCIES_FNAME} \
    {DBS_PARAM} \
    --annotations_fname {MOTIF_ANNOTATIONS_FNAME} \
    --expression_mtx_fname {F_LOOM_SCE} \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 20

#III) AUCELL
!pyscenic aucell \
    {F_LOOM_SCE} \
    reg.csv \
    --output {f_pyscenic_output} \
    --num_workers 20

import json 
import zlib
import base64

# collect SCENIC AUCell output
lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, 'auc2.csv')
auc_mtx.to_csv(AUCELL_MTX_FNAME)

