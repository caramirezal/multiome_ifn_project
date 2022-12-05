## Differential accesibility analysis of human cells pre-treated with (or not) with IFN
## and then stimulated with dsRNA


## Dependencies
library(Seurat)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(viridis)
library(readr)
library(annotatr)
library(GenomicRanges)
library(ggrepel)
library(Matrix)
library(Signac)
library(tibble)


## Initial settings
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')

## Random seed
set.seed(333)


## path to the folder where peaks are stored
path2peaks <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/peaks' 


## path2analysis
path2analysis <- paste0(path2project, 'analysis/differential_accessibility_analysis')
if ( ! dir.exists(path2analysis)){
        dir.create(path2analysis)
}


## path to output figures
path2figures <- paste0(path2project, '/figures/differential_accessibility_analysis')
if ( ! dir.exists(path2figures)){
        dir.create(path2figures)
}



##-----------------------
## Reprocess seurat peaks
reprocess_seurat <- FALSE
if ( reprocess_seurat == TRUE ){
        ##--------------------------
        ## Loading the peak matrix
        peak.mtx <- readMM(file = paste0(path2peaks, '/matrix.mtx'))
        barcodes.tsv <- read_tsv(file = paste0(path2peaks, '/barcodes.tsv'), col_names = FALSE)
        ranges.df <- read_tsv(file = paste0(path2peaks, '/peaks.bed'), FALSE)
        colnames(ranges.df) <- c('seqnames','start', 'end', 'length', 'strand', 'index')
        ranges <- makeGRangesFromDataFrame(ranges.df)
        
        
        
        ##---------------------------
        ## Construction of the annotated reference
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
        broads <- GenomicFeatures::genes(txdb)
        #subject <- broads[ seqnames(broads) %in% seqlevels(ranges) ]
        
        ## Construction of the gene annotations
        TxDb(Homo.sapiens) <- txdb
        tx.list <- transcriptsBy(Homo.sapiens, columns = "SYMBOL")
        tx <- unlist(tx.list)
        tx <- subset(tx, SYMBOL %in% geneSets)
        
        
        
        ##-------------------------
        ## Definition of seurat object
        ## Using seurat for normalization and UMAP projection
        rownames(peak.mtx) <- with(ranges.df, paste0(seqnames, ':', start, '-', end)) 
        colnames(peak.mtx) <- barcodes.tsv$X1
        seurat <- CreateSeuratObject(
                counts = peak.mtx,
                assay = 'ATAC', 
                min.cells = 1, 
                min.features = 1
        )
        seurat$'samples' <- barcodes.tsv$X2 
        
        saveRDS(seurat, 
                file = paste0(path2analysis, 'atac_peaks_seu.rds'), 
                compress = TRUE)
}





## Loading scRNA-Seq data in seurat object
seurat <- read_rds(paste0(path2analysis, 'atac_peaks_seu.rds'))




##------------------------
## Differential peak analysis
sample_comparissons <-  c('4_3h_pIFN_dsRNA', '6_3h_-IFN_dsRNA')
Idents(seurat) <- seurat$samples
markers.dsrna <- FindMarkers(subset(seurat, 
                                    samples %in% sample_comparissons), 
                             ident.1 = sample_comparissons[1], 
                             logfc.threshold = 0.01, 
                             min.pct = 0.01)



##--------------------------
#3 Saving the analysis
if ( ! dir.exists(paths = paste0(path2analysis, '/dea'))){
        dir.create(path = paste0(path2analysis, '/dea'))
}
markers.dsrna %>%
        arrange(desc(avg_log2FC)) %>%
        write_tsv(paste0(path2analysis, 
                         '/dea/diff_peaks_4_3h_pIFN_dsRNA_VS_6_3h_-IFN_dsRNA.tsv.gz'))



##----------------------------
## Getting IFN related genes
## Loading signatures
sign.df <- readxl::read_xlsx('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/signatures/ifn_signatures.xlsx')
## Downloading list of genes
url.list <- list(ifn='https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=HECKER_IFNB1_TARGETS&fileType=txt')
## Loading gene sets
paths.hum.list <- lapply(url.list, readLines)
geneSets <- lapply(paths.hum.list, function(path) path[3:length(path)])
geneSets$'antiviral_response_total' <- sign.df$`antiviral response total`[!is.na(sign.df$`antiviral response total`)]
geneSets$'RIGI_only' <- sign.df$`RIG-I only`[!is.na(sign.df$`RIG-I only`)]
geneSets$'IFN_only' <- sign.df$`IFN only`[!is.na(sign.df$`IFN only`)]
geneSets <- unlist(geneSets) %>% unique()




markers.dsrna %>% 
        ggplot(aes(x=avg_log2FC, y=-log10(p_val))) +
        geom_point() +
        geom_hline(yintercept = -log10(0.05), 
                   linetype='dashed', 
                   colour='red') +
        theme_classic() +
        labs(title = sample_comparissons[1],
             subtitle = paste('Vs ', sample_comparissons[2]))



sel.markers <- markers.dsrna 
dm_regions <- data.frame(seqnames=gsub(':.*', '', rownames(sel.markers)),
                         start=gsub('.*:|-.*', '', rownames(sel.markers)),
                         end=gsub('.*-', '', rownames(sel.markers)),
                         strand='*',
                         p_val=sel.markers$p_val,
                         p_val_adj=sel.markers$p_val_adj, 
                         log2fc=sel.markers$avg_log2FC)
dm_regions <- makeGRangesFromDataFrame(dm_regions, keep.extra.columns = TRUE)

annots = c('hg38_basicgenes')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)

# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
        regions = dm_regions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)
# A GRanges object is returned
print(dm_annotated)


# Coerce to a data.frame
df_dm_annotated = data.frame(dm_annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_dm_annotated))
dim(df_dm_annotated)

df_dm_annotated %>%
        filter(annot.symbol %in% geneSets) %>%
        pull(annot.symbol) %>%
        unique() %>%
        sort()




selected_markers <- arrange(df_dm_annotated, desc(abs(log2fc))) %>%
                        filter(p_val < 0.05) %>%
                        pull(annot.symbol) %>%
                        unique() %>%
                        head(10)
                        


df_dm_annotated <- df_dm_annotated %>% 
        mutate(highlight=ifelse( ( annot.symbol %in% geneSets & p_val < 0.05 ) |
                                        ( annot.symbol == 'ISG15' ), TRUE, FALSE)) %>%
        mutate(gene_label=ifelse(highlight==TRUE, 
                                 annot.symbol,
                                 ''))


df_dm_annotated.unique <- df_dm_annotated
df_dm_annotated.unique <- dplyr::select(df_dm_annotated.unique, seqnames, start, end)
df_dm_annotated.unique <- df_dm_annotated[!duplicated(df_dm_annotated.unique),]
head(df_dm_annotated.unique)

pdf(paste0(path2figures, '/vulcano_plot_treatedVsUntreated_dsRNA.pdf'), 
    height = 12, width = 12)
df_dm_annotated.unique %>%
                ggplot(aes(x=log2fc, y=-log10(p_val),
                           label=gene_label)) +
                geom_point(colour='black') +
                geom_point(data=subset(df_dm_annotated.unique,
                                  annot.symbol %in% geneSets),
                           aes(x=log2fc, y=-log10(p_val),
                               label=gene_label),
                           colour='red') +
                geom_text_repel(max.overlaps = 1000,
                                size=4,
                                force=0.5, 
                                colour='red') +
                geom_hline(yintercept = -log10(0.05), 
                         linetype='dashed', 
                          colour='red') +
                theme_classic() +
                labs(title = sample_comparissons[1],
                     subtitle = paste('Vs ', sample_comparissons[2]),
                     x='Average Log2FC',
                     y='-log10(p_val)')
dev.off()




##----------------------------
## UMAP projection of samples only using features with
## differential peak patterns
sample_comparissons <-  c('4_3h_pIFN_dsRNA', '6_3h_-IFN_dsRNA')
markers.dsrna <- FindMarkers(subset(seurat, 
                                    samples %in% sample_comparissons), 
                             ident.1 = sample_comparissons[1], 
                             logfc.threshold = 0, 
                             min.pct = 0)
selected.markers <- markers.dsrna %>% 
                        filter(p_val < 0.05) %>%
                        rownames()
seurat.sub <- seurat[selected.markers, ]
seurat.sub <- subset(seurat.sub, samples %in% sample_comparissons)
seurat.sub <- RunTFIDF(seurat.sub, assay = 'ATAC')
seurat.sub <- FindTopFeatures(seurat.sub, min.cutoff = "q0")
seurat.sub <- RunSVD(seurat.sub)
seurat.sub <- RunUMAP(seurat.sub, 
                      reduction = 'lsi', 
                      dims=1:20, 
                      n.neighbors = 100, 
                      min.dist = 0.01)
umap.dsRNA <- DimPlot(seurat.sub, 
        group.by = 'samples', 
        cols = sample.colors, 
        pt.size = 3) +
        theme_void()

markers.dsrna %>%
        rownames_to_column('peak') %>%
        write_tsv(file = paste0(path2analysis, 
                                'differential_peaks_dsRNA.tsv.gz'))


##----------------------------
## UMAP projection of samples only using features with
## differential peak patterns
sample_comparissons <-  c('5_3h_pIFN_polyC', '7_3h_-IFN_polyC')
markers.polyC <- FindMarkers(subset(seurat, 
                                    samples %in% sample_comparissons), 
                             ident.1 = sample_comparissons[1], 
                             logfc.threshold = 0, 
                             min.pct = 0)

selected.markers <- markers.polyC %>% 
        filter(p_val < 0.05) %>%
        rownames()
seurat.sub <- seurat[selected.markers, ]
seurat.sub <- subset(seurat.sub, samples %in% sample_comparissons)
seurat.sub <- RunTFIDF(seurat.sub, assay = 'ATAC')
seurat.sub <- FindTopFeatures(seurat.sub, min.cutoff = "q0")
seurat.sub <- RunSVD(seurat.sub)
seurat.sub <- RunUMAP(seurat.sub, 
                      reduction = 'lsi', 
                      dims=1:20, 
                      n.neighbors = 100, 
                      min.dist = 0.01)
umap.polyC <- DimPlot(seurat.sub, 
        group.by = 'samples', 
        cols = sample.colors, 
        pt.size = 3) +
        theme_void()

markers.polyC %>%
        rownames_to_column('peak') %>%
        write_tsv(file = paste0(path2analysis, 
                                'differential_peaks_polyC.tsv.gz'))


pdf(file = paste0(path2figures, '/umap_diff_peaks.pdf'),
    height = 4.5, width = 12)
umap.dsRNA + umap.polyC
dev.off()





##--------------------------
## Checking the order or rows
all(rownames(markers.dsrna) == rownames(markers.polyC))
## Not in the same order, rearranging
markers.dsrna <- markers.dsrna[rownames(markers.polyC),]
## checking again
all(rownames(markers.dsrna) == rownames(markers.polyC))





## Merging
colnames(markers.dsrna) <- paste0(colnames(markers.dsrna), '_dsRNA')
colnames(markers.polyC) <- paste0(colnames(markers.polyC), '_polyC')
markers <- cbind(markers.dsrna, markers.polyC)




## Adding annotations
anns <- read_rds('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/peaks_annotation/annotated_peaks_unambiguously.rds')
rownames(anns) <- with(anns, paste0(seqnames, ':', start, '-', end))
peaks.int <- intersect(rownames(anns), rownames(markers)) 
markers.nd <- markers[!duplicated(peaks.int), ]
anns.nd <- anns[!duplicated(peaks.int), ]
## checking orders
markers.nd <- rownames_to_column(markers.nd, var = 'peak')
anns.nd <- rownames_to_column(anns.nd, var = 'peak')
markers.nd <- merge(markers.nd, anns.nd)
dim(markers.nd)


selected.markers <- filter(markers.nd, p_val_dsRNA<0.05 & p_val_polyC<0.05) %>%
                        pull(annot.symbol)
intersect(selected.markers, geneSets)
                        



##---------------------------
## Highlighting IFN genes with statistical significant change
markers.nd <- markers.nd %>%
        filter(p_val_dsRNA<0.05 | p_val_polyC<0.05) %>%
        mutate(highlight=ifelse(annot.symbol %in% geneSets,
                                TRUE, FALSE)) %>%
        mutate(gene=ifelse(highlight==TRUE, annot.symbol, ''))




pdf(file = paste0(path2figures, '/differential_peak_analysis_IFN(+Vs-)_dsRNa_Vs_polyIC.pdf'),
    height = 8, width = 8)
markers.nd %>%        
        ggplot(aes(x = avg_log2FC_polyC,
                   y = avg_log2FC_dsRNA, 
                   colour = highlight, 
                   label = gene)) +
                geom_point() +
                geom_point(data = subset(markers.nd,
                           annot.symbol %in% degs.polyC),
                           aes(x = avg_log2FC_polyC,
                           y = avg_log2FC_dsRNA),
                           colour='orange') +
                geom_text_repel(max.overlaps = 40,
                                force = 3) +
                geom_hline(yintercept = 0,
                           linetype = 'dashed',
                           colour='red') +
                geom_vline(xintercept = 0,
                           linetype = 'dashed',
                           colour='red') +
                theme_classic() +
                theme(legend.position = 'none') +
                scale_color_manual(values = c('black', 'red')) +
                labs(x='Average Log2FC IFN(+/-) | polyC',
                     y='Average Log2FC IFN(+/-) | dsRNA') 
dev.off()

