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
# python libs
import decoupler as dc
import session_info

# LIANA loading
import liana as li
from liana.mt import rank_aggregate
?rank_aggregate.__call__
rank_aggregate.describe()

# import all individual methods
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean


# Cells at recurrence
adata_recurrent = adata[adata.obs['Stage'] == 'Recurrent'].copy()

# run cellphonedb
cellphonedb(adata_recurrent,
            groupby='custom_annotation', 
            resource_name='consensus',
            expr_prop=0.1,
            verbose=True, key_added='cpdb_res')

adata_recurrent.uns['cpdb_res'].head()
adata_recurrent.uns['cpdb_res'].to_csv(os.path.join(datadir,'cpdb_res.csv'))

li.pl.dotplot(adata = adata_recurrent, 
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True, # we inverse sign since we want small p-values to have large sizes
              source_labels=['Oligodendrocytes'],
              target_labels=['Neurons','Immune cells'],
              figure_size=(8, 7),
              # finally, since cpdbv2 suggests using a filter to FPs
              # we filter the pvals column to <= 0.05
              filter_fun=lambda x: x['cellphone_pvals'] <= 0.05,
              uns_key='cpdb_res' # uns_key to use, default is 'liana_res' 
             )

li.plotting.tileplot(adata = adata_recurrent, 
                         fill='means',
                         label='props',
                         label_fun=lambda x: f'{x:.2f}',
                         top_n=25, 
                         orderby='cellphone_pvals',
                         orderby_ascending=True,
                         source_labels='Oligodendrocytes',
                         target_labels=['Neurons','Immune cells'],
                         uns_key='cpdb_res', # NOTE: default is 'liana_res'
                         source_title='Ligand',
                         target_title='Receptor',
                         figure_size=(8, 7)
                         )

li.plotting.tileplot(adata = adata_recurrent, 
                         fill='means',
                         label='props',
                         label_fun=lambda x: f'{x:.2f}',
                         top_n=25, 
                         orderby='cellphone_pvals',
                         orderby_ascending=True,
                         source_labels='Immune cells',
                         target_labels=['Tumor','Oligodendrocytes','Neurons'],
                         uns_key='cpdb_res', # NOTE: default is 'liana_res'
                         source_title='Ligand',
                         target_title='Receptor',
                         figure_size=(8, 7)
                         )

# Run rank_aggregate
li.mt.rank_aggregate(adata_recurrent, 
                     groupby='custom_annotation',
                     resource_name='consensus',
                     expr_prop=0.1,
                     verbose=True)

li.pl.circle_plot(adata_recurrent,
                  groupby='custom_annotation',
                  score_key='magnitude_rank',
                  inverse_score=True,
                  source_labels='Oligodendrocytes',
                  filter_fun=lambda x: x['specificity_rank'] <= 0.05,
                  pivot_mode='counts', # NOTE: this will simply count the interactions, 'mean' is also available
                  figure_size=(10, 10),
                  )

my_plot = li.pl.dotplot(adata = adata_recurrent, 
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
                        source_labels=['Neurons', 'Oligodendrocytes', 'Immune cells','Tumor'],
                        target_labels=['Neurons', 'Oligodendrocytes','Immune cells'],
                        filter_fun=lambda x: x['specificity_rank'] <= 0.01
                       )
my_plot

li.pl.circle_plot(adata_recurrent,
                  groupby='custom_annotation',
                  score_key='magnitude_rank',
                  inverse_score=True,
                  source_labels=['Oligodendrocytes','Neurons','Immune cells'],
                  filter_fun=lambda x: x['specificity_rank'] <= 0.05,
                  pivot_mode='counts', # NOTE: this will simply count the interactions, 'mean' is also available
                  figure_size=(10, 10),
                  )

# we import plotnine
import plotnine as p9

(my_plot +
 # change theme
 p9.theme_dark() +
 # modify theme
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11),
     axis_text_x=p9.element_text(rotation=45, ha='right'),
     figure_size=(12, 15)
 )
)

li.pl.dotplot(adata = adata_recurrent, 
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels='Oligodendrocytes',
              target_labels=['Neurons','Immune cells'],
              top_n=25, 
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(8, 7)
             )

# Cells at primary condition

adata_primary = adata[adata.obs['Stage'] == 'Primary'].copy()
# run cellphonedb
cellphonedb(adata_primary,
            groupby='custom_annotation', 
            resource_name='consensus',
            expr_prop=0.1,
            verbose=True, key_added='cpdb_prm')

li.plotting.tileplot(adata = adata_primary, 
                         fill='means',
                         label='props',
                         label_fun=lambda x: f'{x:.2f}',
                         top_n=25, 
                         orderby='cellphone_pvals',
                         orderby_ascending=True,
                         source_labels='Oligodendrocytes',
                         target_labels=['Neurons','Immune cells','Tumor'],
                         uns_key='cpdb_prm', # NOTE: default is 'liana_res'
                         source_title='Ligand',
                         target_title='Receptor',
                         figure_size=(8, 7)
                         )

adata_primary.uns['cpdb_prm'].to_csv(os.path.join(datadir,'cpdb_primary.csv'))

li.plotting.tileplot(adata = adata_primary, 
                         fill='means',
                         label='props',
                         label_fun=lambda x: f'{x:.2f}',
                         top_n=25, 
                         orderby='cellphone_pvals',
                         orderby_ascending=True,
                         source_labels='Oligodendrocytes',
                         target_labels=['Neurons','Immune cells'],
                         uns_key='cpdb_prm', # NOTE: default is 'liana_res'
                         source_title='Ligand',
                         target_title='Receptor',
                         figure_size=(8, 7)
                         )

li.mt.rank_aggregate(adata_primary, 
                     groupby='custom_annotation',
                     resource_name='consensus',
                     expr_prop=0.1,
                     verbose=True)

primary=li.pl.dotplot(adata = adata_primary, 
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              source_labels=['Oligodendrocytes', 'Neurons', 'Immune cells','Tumor'],
              target_labels=['Oligodendrocytes', 'Neurons','Immune cells'],
              orderby='magnitude_rank',
              filter_fun=lambda x: x['specificity_rank'] <= 0.01,
              orderby_ascending=True,
              figure_size=(10, 12)
             )

(primary +
 # change theme
 p9.theme_dark() +
 # modify theme
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11),
     axis_text_x=p9.element_text(rotation=45, ha='right'),
     figure_size=(12, 15)
 )
)

li.pl.circle_plot(adata_primary,
                  groupby='custom_annotation',
                  score_key='magnitude_rank',
                  inverse_score=True,
                  source_labels='Oligodendrocytes',
                  filter_fun=lambda x: x['specificity_rank'] <= 0.05,
                  pivot_mode='counts', # NOTE: this will simply count the interactions, 'mean' is also available
                  figure_size=(10, 10),
                  )