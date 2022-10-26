## This script contains the analysis of the IFN signatures
## using multiomic data

## Dependencies
library(Seurat)
library(viridis)
library(AUCell)

## Initial settings
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')


## Loading pre-processed seurat object
## the object was generated using the script stored in:
## scripts/rna_seq_processed_individually.R 
seurat <- read_rds(paste0(path2project, 'data/seurat/ifn_treated_cells_rnaseq_processed_seurat.rds'))



## Downloading list of genes
url.list <- list(il12_Down='https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=GSE15930_NAIVE_VS_24H_IN_VITRO_STIM_IL12_CD8_TCELL_DN&fileType=txt',
                 il12_Up='https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=GSE15930_NAIVE_VS_24H_IN_VITRO_STIM_IL12_CD8_TCELL_UP&fileType=txt',
                 ifn='https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=HECKER_IFNB1_TARGETS&fileType=txt',
                 TCR='https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY&fileType=txt')


## Loading gene sets
paths.hum.list <- lapply(url.list, readLines)
geneSets <- lapply(paths.hum.list, function(path) path[3:length(path)])


lapply(geneSets, head)

## Calculating the scores
cells_rankings <- AUCell_buildRankings(seurat@assays$SCT@scale.data)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.05)
auc.mtx <- getAUC(cells_AUC)
auc.df <- as.data.frame(t(auc.mtx))
head(auc.df)


##-----------------
## Adding information to the seurat object
all(colnames(seurat) == rownames(auc.df) )
seurat$'ifn_aucell_signature' <- auc.df$ifn


FeaturePlot(seurat, 
            features = 'ifn_aucell_signature', 
            pt.size = 3) +
                scale_color_viridis()
