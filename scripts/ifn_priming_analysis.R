## This script contains the analysis of the primed state after IFN pre-treatment

## dependencies
library(dplyr)
library(readr)
library(ggplot2)
library(Seurat)
library(ggrepel)
library(ggvenn)
library(UpSetR)


## Loading initial settings
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')



## Load differential peaks between IFN treated Vs non-treated in the
## non polyC context
dif.peaks <- read_tsv('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/differential_accessibility_analysis/differential_peaks_polyC.tsv.gz')
head(dif.peaks)



## Differential peak analysis
seurat <- readRDS('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/differential_accessibility_analysis/atac_peaks_seu.rds')


##------------------------
## Differential peak analysis
sample_comparissons <-  c('5_3h_pIFN_polyC', '7_3h_-IFN_polyC')
Idents(seurat) <- seurat$samples
markers <- FindMarkers(subset(seurat, 
                              samples %in% sample_comparissons), 
                       ident.1 = sample_comparissons[1], 
                       logfc.threshold = 0, 
                       min.pct = 0)
dim(markers)


## Annotation of peaks
peak.anns <- read_rds('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/peaks_annotation/annotated_peaks_unambiguously.rds') 
peak.anns <- mutate(peak.anns, peak=paste0(seqnames, ':', start, '-', end))
markers <- merge(markers, peak.anns)


annot.type.df <- peak.anns$annot.type %>%
        table() %>%
        as.data.frame()
colnames(annot.type.df) <- c('region_type', 'frequency')
barplot.anns <- annot.type.df %>%
        ggplot(aes(x=region_type,y=frequency)) +
                geom_bar(position = 'stack',
                         stat = 'identity',
                         fill='steelblue') +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                                         size = 10),
                              axis.text.y = element_text(size = 10)) +
                        labs(x='', y='Frequency')

table(peak.anns$annot.type) %>% sort()
intersect(colnames(peak.anns), colnames(markers))
dim(markers)
head(markers)
sum(is.na(markers$annot.symbol))



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






## Vulcano plot
markers$'peak' <- rownames(markers)
markers <- markers %>%
        arrange(desc(avg_log2FC)) %>%
        mutate(ranking=1:nrow(markers)) %>%
        mutate(highlight=ifelse(annot.symbol %in% geneSets,
                                TRUE, FALSE)) %>%
        mutate(gene_label=ifelse( ( annot.symbol == 'ISG15' ) | 
                                          ( annot.symbol %in% geneSets &
                                         p_val < 0.05 ),
                                 annot.symbol, ''))
markers %>%
        ggplot(aes(x=avg_log2FC, y=-log10(p_val),
                   label=gene_label)) +
                geom_point() +
                geom_point(data=subset(markers, highlight==TRUE), 
                           aes(x=avg_log2FC, y=-log10(p_val)),
                           colour='red') +
                geom_hline(yintercept = -log10(0.05),
                           linetype='dashed', 
                           colour='red') +
                geom_text_repel(max.overlaps = 100, 
                                colour='red', 
                                force = 1) +
                theme_classic() +
                labs(x='Average Log2 FC',
                     y='-Log10(P-Value)')



vulcano.dif.peaks <- markers %>%
        ggplot(aes(x=avg_log2FC, y=-log10(p_val))) +
        geom_point() +
        geom_point(data=subset(markers, p_val < 0.05 ), 
                   aes(x=avg_log2FC, y=-log10(p_val)),
                   colour='red') +
        geom_hline(yintercept = -log10(0.05),
                   linetype='dashed', 
                   colour='red') +
        theme_classic() +
        labs(x='Average Log2 FC',
             y='-Log10(P-Value)', 
             title = 'Gained peaks in IFN pre- Vs Non treated cells',
             subtitle = 'PolyC')
        


gained.peaks <- subset(markers, p_val < 0.05 & 0 < avg_log2FC ) %>% pull(`annot.symbol`)
losed.peaks <- subset(markers, p_val < 0.05 & avg_log2FC < 0 ) %>% pull(`annot.symbol`)


## Loading differentially expressed genes comparing dsRNA Vs polyC
## for the two conditions pre-treated and non-treated
degs <- read.table('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/ifn_signatures/dsRNAVspolyC_across_conditions.tsv')
str(degs)

degs.up.ifn.pre <- filter(degs, p_val_pIFN < 0.05 & avg_log2FC_pIFN > 0 ) %>% pull(gene)
degs.up.ifn.non <- filter(degs, p_val_noIFN < 0.05 & avg_log2FC_noIFN > 0 ) %>% pull(gene)
degs.down.ifn.pre <- filter(degs, p_val_pIFN < 0.05 & avg_log2FC_pIFN < 0 ) %>% pull(gene)
degs.down.ifn.non <- filter(degs, p_val_noIFN < 0.05 & avg_log2FC_noIFN < 0 ) %>% pull(gene)


gene_list <- list(gained_peaks=gained.peaks,
                  losed_peaks=losed.peaks,
                  degs_up_pre_IFN=degs.up.ifn.pre,
                  degs_up_no_IFN=degs.up.ifn.non,
                  degs_down_pre_IFN=degs.down.ifn.pre,
                  degs_down_no_IFN=degs.down.ifn.non)
gene_list <- lapply(gene_list, unique)


upset_plot <- upset(fromList(gene_list), nsets = 100)

ggvenn(data = list(degs_up_pre_IFN=degs.up.ifn.pre,
                   degs_up_no_IFN=degs.up.ifn.non),
  )


barplot.anns



degs <- degs %>% mutate(regulation_type='Non-regulated')
degs$regulation_type[degs$gene %in% degs.up.ifn.pre] <- 'Up-regulated in pre-treated'
degs$regulation_type[degs$gene %in% degs.up.ifn.non] <- 'Up-regulated in non-treated'
degs$regulation_type[degs$gene %in% intersect(degs.up.ifn.non, degs.up.ifn.pre)] <- 'Up-regulated shared'
degs <- degs %>% mutate(increased_peak='No-increase')
degs$increased_peak[degs$gene %in% gained.peaks$gene_label] <- 'Increased peak'

table(degs$increased_peak, degs$regulation_type)
subset(degs, regulation_type == 'Up-regulated in pre-treated' 
                & increased_peak == 'Increased peak' )
