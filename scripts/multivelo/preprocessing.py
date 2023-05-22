## This script contains the analysis of human cell lines primed with interferon
## and treated or not dsRNA


##  Dependencies
import collections 
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import loompy
import matplotlib
import anndata


## Initial settings
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)


## defining folder to store figures
path2figures = '/media/ag-cherrmann/cramirez/multiome_ifn_project/figures/multivelo'
if not os.path.exists(path2figures):
    os.mkdir(path2figures)



##------------------------------------------------------------------------------------------
## Reading in unspliced and spliced counts




## If there more than 1 Loom files they need to be merged
path2loom = '/media/ag-cherrmann/cramirez/multiome_ifn_project/data/multivelo/loom/'
loom_files = os.listdir(path2loom)
abs_paths = [ path2loom + file for file in loom_files]
combined_loom_path = path2loom + '/multiome_ifn_combined.loom'
if not os.path.exists(combined_loom_path):
    loompy.combine(
    abs_paths, 
    combined_loom_path, 
    key='Accession')
    ## cleaning loom intermediate files
    for file in abs_paths: os.remove(file)



## reading combined data
adata = sc.read_loom(combined_loom_path)


## Checking adata content
adata
adata.obs_names
adata.var_names
adata.var



##------------------------------------------------------------------------------------------
## Velocyto pipeline


## Removing condition label from barcode
adata.obs['barcode']  = [ st[(st.index(':')+1):] for st in adata.obs_names ] 
adata.obs


## Loading annotations
anns_path = '/media/ag-cherrmann/cramirez/multiome_ifn_project/data/seurat/seurat_annotations.tsv.gz'
anns = pd.read_csv(anns_path, sep='\t')
anns['barcode'] = [ st[ :(st.index('-'))] for st in anns['barcode'] ] 
anns.barcode



## Intersecting barcodes with preprocessed and annotated cells
barcodes_adata = adata.obs['barcode'].tolist()
barcodes_anns = anns['barcode'].tolist()


## Removing duplicated barcodes
non_dup_adata_bc = [i for i in barcodes_adata if barcodes_adata.count(i) == 1 ]
non_dup_adata_anns = [i for i in barcodes_anns if barcodes_anns.count(i) == 1 ]
barcodes_intersection = set(non_dup_adata_bc).intersection(set(non_dup_adata_anns))
barcodes_intersection = list(barcodes_intersection)



## Selecting cells in adata
adata = adata[adata.obs['barcode'].isin(barcodes_intersection),]
anns_processed = anns[anns['barcode'].isin(barcodes_intersection)]



## Rearranging annotation dataframe
reordering = pd.Series(anns_processed['barcode'], dtype='category')
reordering = reordering.cat.reorder_categories(adata.obs['barcode'], ordered=True)
anns_processed.barcode = reordering
anns_processed.sort_values(by='barcode', inplace=True)



## Checkin order in barcodes
len(anns_processed['barcode']) == len(adata.obs['barcode'])
matching_barcode = [ anns_processed['barcode'].iloc[i] == adata.obs['barcode'].iloc[i]  for i in range(len(anns_processed['barcode']))]
sum(matching_barcode) == len(matching_barcode)


## Adding annotations to the adata object
anns_processed = anns_processed.set_axis(anns_processed.barcode, axis='index')
adata.obs = anns_processed
## Adding UMAP information to the adata object
adata.obsm['X_umap'] = adata.obs[['UMAP_1', 'UMAP_2']].to_numpy()
adata.obs.columns


## Normalization
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)



## Pseudotime calculation
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)


## Plotting
color_mapping = {'7_3h_-IFN_polyC' : 'lime',
                   '6_3h_-IFN_dsRNA' : 'limegreen',
                   '3_16h_-IFN_dsRNA' : 'darkgreen',
                   '5_3h_pIFN_polyC' : 'lightsalmon',
                   '4_3h_pIFN_dsRNA' : 'orange'}

figure_file_name = path2figures + '/velocity_umap.pdf'
## Plot settings
matplotlib.rcParams['figure.figsize'] = 6, 6
scv.pl.velocity_embedding_stream(adata, 
                                 basis='umap', 
                                 color = 'orig.ident', 
                                 palette=color_mapping)
plt.axis('off')
plt.savefig(figure_file_name , 
            format='pdf',
            bbox_inches='tight', 
            pad_inches=0.25)








##-----------------------------------------
## ATAC processing



## Loading ATAC seq peak matrix
path2peaks = '/media/ag-cherrmann/cramirez/multiome_ifn_project/data/peaks/'
counts = scipy.io.mmread(path2peaks + 'matrix.mtx').T
#counts = scipy.sparse.csr_matrix(counts)
regions = pd.read_csv(path2peaks + 'peaks.bed', sep='\t', header=None, names=['chr','start','end', 'len', 'dir', 'index'])
cells = pd.read_csv(path2peaks + 'barcodes.tsv', sep='\t', header=None, names=['barcodes', 'condition'])

# then initialize a new AnnData object
adata_atac = anndata.AnnData(X=counts, obs=cells, var=regions)


## Aggregate peaks
adata_atac = mv.aggregate_peaks_10x(adata_atac,
                                    '/media/ag-cherrmann/cramirez/multiome_ifn_project/analysis/peaks_annotation/peak_annotation.tsv',
                                    '/media/ag-cherrmann/cramirez/multiome_ifn_project/data/multivelo/feature_linkage.bedpe',
                                    verbose=True)


 



## visualization of count distribution
counts_file_name = path2figures + '/counts_quality_control.pdf'
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
plt.style.use('default')
#matplotlib.rcParams['figure.figsize'] = 6, 6
plt.hist(adata_atac.X.sum(1), bins=30, range=(0, 2));
plt.savefig(counts_file_name , 
            format='pdf')
            
