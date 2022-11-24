## This script contains the analysis of ATAC-Seq of human cells pre-treated or not
## with IFN and stimulated with polyC

## Dependencies
library(Seurat)
library(Matrix)
library(dplyr)
library(readr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Homo.sapiens)
library(IRanges)
library(GenomicRanges)
library(Signac)
library(cluster)


## Setting initial parameters
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')



## path to store output figures
path2figures <- paste0(path2project, '/figures/ifn_focused_analysis.R')


## path to the folder where peaks are stored
path2peaks <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/peaks' 



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




##--------------------------
## Loading the peak matrix
peak.mtx <- readMM(file = paste0(path2peaks, '/matrix.mtx'))
barcodes.tsv <- read_tsv(file = paste0(path2peaks, '/barcodes.tsv'), col_names = FALSE)
ranges.df <- read_tsv(file = paste0(path2peaks, '/peaks.bed'), FALSE)
colnames(ranges.df) <- c('seqnames','start', 'end', 'length', 'strand', 'index')
ranges <- makeGRangesFromDataFrame(ranges.df)




## Writing peaks info
write_tsv(ranges.df, 
          file = paste0(path2peaks, '/peaks.bed'), 
          col_names = FALSE)





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
## Getting the regions that overlaps (inside gene regions)
## peaks_in_genes
overlaps.gene <- findOverlaps(query = ranges, subject = tx)
peaksInGenes.indxs <- queryHits(overlaps.gene) %>% unique()
length(peaksInGenes.indxs)




##-------------------------
## Regions around 10 kb of the TSS
#tss.peaks <- resize(tx, width = 1, fix = 'start')
tss.peaks <- tx
window.width <- 100000 
start(tss.peaks) <- start(tss.peaks) - window.width
end(tss.peaks) <- end(tss.peaks) + window.width
overlaps.gene <- findOverlaps(query = ranges, subject = tss.peaks)
peaksInTSS.indxs <- queryHits(overlaps.gene) %>% unique()
length(peaksInTSS.indxs)


## List of peaks subsets
peaks.subsets.list <- list('peaks_in_genes'=peaksInGenes.indxs,
                           'peaks_tss'=peaksInTSS.indxs)




##-------------------------
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




##------------------------------
## Subsetting the peak matrix
#DefaultAssay(seurat) <- 'ATAC'
seurat <- seurat[peaks.subsets.list$peaks_tss, ]
seurat <- subset(seurat, samples %in% c('4_3h_pIFN_dsRNA', '3_16h_-IFN_dsRNA'))
seurat <- RunTFIDF(seurat, assay = 'ATAC')
seurat <- FindTopFeatures(seurat, min.cutoff = "q0")
seurat <- RunSVD(seurat)
seurat <- RunUMAP(seurat, reduction = 'lsi', dims=1:50)


DimPlot(seurat, reduction = 'umap', group.by = 'samples', split.by = 'samples',
        cols = sample.colors)




##-------------------------
## Calculation of Silhouette to address separability of clusters
head(seurat@reductions$umap@cell.embeddings)
dists <- dist(seurat@reductions$umap@cell.embeddings)
sil <- silhouette(x = as.integer(factor(seurat$samples)), 
           dist = dists)
mean(sil[, 3])

