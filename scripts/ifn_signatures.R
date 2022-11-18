## This script contains the analysis of the IFN signatures
## using multiomic data

## Dependencies
library(Seurat)
library(viridis)
library(AUCell)
#library(tidyverse)
library(ggrepel)
library(dplyr)
library(tibble)

set.seed(333)

dev.off()

## Initial settings
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')


## path to output figures
path2figures <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/ifn_signatures'
if ( ! dir.exists(path2figures)){
        dir.create(path2figures)
}

## Loading pre-processed seurat object
## the object was generated using the script stored in:
## scripts/rna_seq_processed_individually.R 
seurat <- readRDS(paste0(path2project, 'data/seurat/ifn_treated_cells_rnaseq_processed_seurat.rds'))



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



##----------------
## Visualization

## By sample
sample.colors <- sample.colors[! names(sample.colors) %in% c('1_16h_pIFN_dsRNA','2_16h_pIFN_polyC')]
treatment.colors <- sample.colors
names(treatment.colors) <- sapply(names(sample.colors), function(x) substr(x, 3, nchar(x)))
dev.off()
sample.plot <- DimPlot(seurat, reduction = "umap", 
                       group.by = 'treatment', 
                       cols = treatment.colors, 
                       label = TRUE, 
                       repel = TRUE, 
                       pt.size = 2) + 
                        theme_void() + 
                        labs(title = 'By Sample')

umap.ifn.sign <- FeaturePlot(seurat, 
            features = 'ifn_aucell_signature', 
            pt.size = 3) +
                scale_color_viridis() +
                theme_void() +
                labs(title='IFN Signature')


sample.plot + umap.ifn.sign




##----------------------
## Differential expression analysis by condition
Idents(seurat) <- seurat$orig.ident

recalculate_degs <- FALSE
if ( recalculate_degs == TRUE ) {
        degs <- FindAllMarkers(seurat, 
                               logfc.threshold = 0, 
                               min.pct = 0, 
                               min.cells.feature = 1)
        degs.conditions.list <- split(degs, f = degs$cluster)
        degs.conditions.list <- lapply(degs.conditions.list, 
                                       function(degs) { arrange(degs, 
                                                                desc(avg_log2FC))
                                       })
        path2degs <- paste0(path2project, 'data/degs')
        
        if ( ! dir.exists(path2degs) ){
                dir.create(path2degs) 
                writexl::write_xlsx(
                        degs.conditions.list, 
                        path = paste0(path2degs, '/degs_all_conditions_vs_rest.xlsx')
                ) 
        }
}




################################################################################


##-----------------------
## DEA comparing synergystic genes vs dsRNA upregulated genes
## Comparing 3h IFN+polyIC vs IFN-polyIC
degs.pifn <- FindMarkers(object = subset(seurat, 
                                     orig.ident %in%  c('4_3h_pIFN_dsRNA',
                                                        '5_3h_pIFN_polyC') ), 
                             ident.1 = '4_3h_pIFN_dsRNA', 
                             logfc.threshold = 0,
                             min.pct = 0)
colnames(degs.pifn) <- paste0(colnames(degs.pifn), '_+IFN')
degs.pifn <- rownames_to_column(degs.pifn, var = 'gene')
degs.mifn <- FindMarkers(object = subset(seurat, 
                                          orig.ident %in%  c('6_3h_-IFN_dsRNA',
                                                             '7_3h_-IFN_polyC') ), 
                          ident.1 = '6_3h_-IFN_dsRNA', 
                          logfc.threshold = 0,
                          min.pct = 0)
colnames(degs.mifn) <- paste0(colnames(degs.mifn), '_-IFN')
degs.mifn <- rownames_to_column(degs.mifn, var = 'gene')



##----------------------
## checking up-regulation direction
up.genes.pifn <- pull(head(degs.pifn, 10), var = 'gene')
DotPlot(subset(seurat, orig.ident %in%  c('5_3h_pIFN_polyC', 
                                          '7_3h_-IFN_polyC')), 
        features = up.genes.pifn)

## 
degs.merged <- merge(degs.pifn, degs.mifn)

pdf(paste0(path2figures, '/scat_plot_comparing_IFN_signatures.pdf'),
    height = 6, width = 8)
degs.merged %>%
        mutate(gene_label=ifelse(gene %in% geneSets$ifn, gene, '')) %>%
        ggplot(aes(x=`avg_log2FC_-IFN`,
                   y=`avg_log2FC_+IFN`,
                   size=-log10(`p_val_-IFN`*`p_val_+IFN`), 
                   colour=`pct.1_-IFN`*`pct.1_+IFN`,
                   label=gene_label)) +
                geom_point() +
                geom_text_repel(size= 3, color='black', force = 1) +
                geom_abline(slope = 1, intercept = c(0,0), linetype ='dashed') + 
                geom_hline(yintercept = 0, linetype = 'dashed') +
                geom_vline(xintercept = 0, linetype = 'dashed') +
                scale_color_viridis() +
                theme_classic() +
                theme(legend.position = 'bottom') +
                labs(x='Log2FC dsRNA(+/-) | IFN(-)',
                     y='Log2FC dsRNA(+/-) | IFN(+)')
dev.off()


###########################################################################

##-----------------------
## DEA comparing synergystic genes vs dsRNA upregulated genes
## Comparing 3h IFN+polyIC vs IFN-polyIC
degs.ifnpolyC <- FindMarkers(object = subset(seurat, 
                                             orig.ident %in%  c('5_3h_pIFN_polyC', 
                                                                '7_3h_-IFN_polyC') ), 
                             ident.1 = '5_3h_pIFN_polyC', 
                             logfc.threshold = 0,
                             min.pct = 0)
colnames(degs.ifnpolyC) <- paste0(colnames(degs.ifnpolyC), '_ifnpolyC')
degs.ifnpolyC <- rownames_to_column(degs.ifnpolyC, var = 'gene')
degs.ifndsrna <- FindMarkers(object = subset(seurat, 
                                             orig.ident %in%  c('4_3h_pIFN_dsRNA', 
                                                                '6_3h_-IFN_dsRNA') ), 
                             ident.1 = '4_3h_pIFN_dsRNA', 
                             logfc.threshold = 0,
                             min.pct = 0)
colnames(degs.ifndsrna) <- paste0(colnames(degs.ifndsrna), '_ifndsrna')
degs.ifndsrna <- rownames_to_column(degs.ifndsrna, var = 'gene')



##----------------------
## checking up-regulation direction
up.genes.ifnpolyC <- pull(head(degs.ifnpolyC, 10), var = 'gene')
DotPlot(subset(seurat, orig.ident %in%  c('5_3h_pIFN_polyC', 
                                          '7_3h_-IFN_polyC')), 
        features = up.genes.ifnpolyC)

## 
degs.merged <- merge(degs.ifnpolyC, degs.ifndsrna)


pdf(paste0(path2figures, '/scat_plot_comparing_IFN_signatures_preVsNo_treatment.pdf'),
    height = 6, width = 8)
degs.merged %>%
        mutate(gene_label=ifelse(gene %in% geneSets$ifn, gene, '')) %>%
        ggplot(aes(x=avg_log2FC_ifnpolyC,
                   y=avg_log2FC_ifndsrna,
                   size=-log10(p_val_ifnpolyC*p_val_ifndsrna), 
                   colour=pct.1_ifnpolyC*pct.1_ifndsrna,
                   label=gene_label)) +
        geom_point() +
        geom_text_repel(size= 3, color='black', force = 1) +
        geom_abline(slope = 1, intercept = c(0,0), linetype ='dashed') + 
        geom_hline(yintercept = 0, linetype = 'dashed') +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        scale_color_viridis() +
        theme_classic() +
        theme(legend.position = 'bottom') +
        labs(x='Log2FC IFN (+/-) | polyC',
             y='Log2FC IFN (+/-) | dsRNA')
dev.off()





###########################################################################
## GEX visualization


pdf(file = paste0(path2figures, '/ifn_signatures_by_sample.pdf'),
    width = 20,
    height = 4)
ifn_genes <- geneSets$ifn[geneSets$ifn %in% rownames(seurat)]
DotPlot(seurat, 
        group.by = 'treatment', 
        features = ifn_genes, 
        cols = 'RdYlGn') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x='', y='')
dev.off()


