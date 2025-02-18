# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # Import of matrix data, filtering and Quality Control using Scanpy following the Parse Bioscienes tutorial
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loading libraries and setting the location path for analysis data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

import numpy as np
import pandas as pd
import scanpy as sc
#import scipy
#import os
#import scipy.io as sio
import argparse


# Argumanets
parser = argparse.ArgumentParser()

parser.add_argument("directory", action="store", type=str, metavar='directory', help='Path to the split-pipe output directory that needs processing')
parser.add_argument("outdir", action="store", type=str, metavar='outdir', help='Write result to here')
parser.add_argument("--pipeline", action="store", type=str, metavar='pipeline', help='Script being used as part of pipeline - if so use custom output directory location')

args = parser.parse_known_args()    #Use parse_known_arg to differentiate between arguments pre-specified and those that are not

options = args[0]
data_path = options.directory + '/'
#data_path = 'result_copy/splitpipe_combined_sublibraries/SampleA/DGE_unfiltered/'

# oras://community.wave.seqera.io/library/pip_igraph_leidenalg_scanpy:ae205e6e609e6292

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Adjusting Scanpy default settings
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

sc.settings.verbosity = 1 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=300, figsize=(5,4), format='png')
#sc.settings.figdir = '/mnt/output_figures/'


if(options.pipeline is not None):
    elements = options.outdir.split('/')
    sc.settings.figdir = f'scanpy_singlecellqc_{options.pipeline}_' + '/'.join(elements[:-1])
else:
    sc.settings.figdir = options.outdir


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Reading in data
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# The DGE_filtered folder contains the expression matrix, genes, and files
# NOTE: split-pipe versions older than 1.1.0 used 'DGE.mtx'
#print(data_path)

try:
    adata = sc.read_mtx(data_path + 'count_matrix.mtx')

    # reading in gene and cell data
    gene_data = pd.read_csv(data_path + 'all_genes.csv')
    cell_meta = pd.read_csv(data_path + 'cell_metadata.csv')
except:
  print("Could not import data!")
  exit(0)


try:
    # find genes with nan values and filter
    gene_data = gene_data[gene_data.gene_name.notnull()]
    notNa = gene_data.index
    notNa = notNa.to_list()

    # remove genes with nan values and assign gene names
    adata = adata[:,notNa]
    adata.var = gene_data
    adata.var.set_index('gene_name', inplace=True)
    adata.var.index.name = None
    adata.var_names_make_unique()

    # add cell meta data to anndata object
    adata.obs = cell_meta
    adata.obs.set_index('bc_wells', inplace=True)
    adata.obs.index.name = None
    adata.obs_names_make_unique()

    # Returns the dimensions of the expression matrix (cells, genes)
    print(adata.shape)
except:
      print("Could not filter data!")
      exit(0)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Cell quality control
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

try:
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Scanpy will prepend the string in the save argument with "violin"
    # and save it to our figure directory defined in the first step.
    sc.pl.violin(adata, ['n_genes_by_counts'], save='_n_genes', jitter=0.4)
    sc.pl.violin(adata, ['total_counts'], save='_total_counts', jitter=0.4)
    sc.pl.violin(adata, ['pct_counts_mt'], save='_mito_pct', jitter=0.4)

    # Filter the data
    adata = adata[adata.obs.n_genes_by_counts < 5000,:]
    adata = adata[adata.obs.total_counts < 20000,:]
    adata = adata[adata.obs.pct_counts_mt < 15,:]
    adata.shape # Checking number of cells remaining

    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_gene_vs_transcript_counts')
    print('median transcript count per cell: ' + str(adata.obs['tscp_count'].median(0)))
    print('median gene count per cell: ' + str(adata.obs['gene_count'].median(0)))

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
    sc.pl.highly_variable_genes(adata, save='') # scanpy generates the filename automatically

    # Save raw expression values before variable gene subset
    adata.raw = adata

    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts'])
    sc.pp.scale(adata, max_value=10)
except:
    print("Could not perform cell quality control")
    exit(0)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Principal component analysis
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

try:
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50, save='') # scanpy generates the filename automatically


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # UMAP and Leiden Clustering
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save='_leiden')


    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    # Finding cluster markers
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

    # The head function returns the top n genes per cluster
    top_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
    print(top_markers)
except:
    print("Could not perform dimensionality reduction")


print('Done')