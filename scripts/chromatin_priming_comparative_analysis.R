## This script contains a comparative analysis of peaks gained after IFN treatment
## and the ones gained after dsRNA stimulation

## dependencies
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(readr)
library(viridis)
library(ggrepel)

source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')


## colors
priming_categories <- c('primed'='darkorange',
                        'non-primed'='darkolivegreen3')

rerun_analysis <- FALSE


## Loading seurat object containing peaks
seurat <- readRDS('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/differential_accessibility_analysis/atac_peaks_seu.rds')
## seurat[[]]$samples %>%
##         table()


path2analysis <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/chromatin_priming_comparative_analysis'
if ( ! dir.exists(path2analysis)) {
        dir.create(path2analysis)
}

path2figures <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/chromatin_priming_comparative_analysis'
if ( ! dir.exists(path2figures) ){
        dir.create(path2figures)
}

## DEGs comparing 3 hrs IFN vs non-IFN (during polyC)
Idents(seurat) <- seurat$samples

if ( rerun_analysis == TRUE ) {
        degs.priming <- FindMarkers(
                object = subset(seurat,
                                samples %in% c('7_3h_-IFN_polyC', '5_3h_pIFN_polyC')) ,
                ident.1 = '5_3h_pIFN_polyC',
                ident.2 = '7_3h_-IFN_polyC',
                min.pct = 0,
                logfc.threshold = 0,
                test.use = "wilcox"
        ) 
        colnames(degs.priming) <- paste0(colnames(degs.priming), '_pVsnIFN_polyC')
        degs.priming <- rownames_to_column(degs.priming, var = 'peak')
        head(degs.priming)
        
        write_tsv(degs.priming, 
                  file = paste0(path2analysis,
                                '/diff_peaks_7_3h_-IFN_polyC_VS_5_3h_pIFN_polyC.tsv.gz'))
}



degs.priming <- read_tsv(paste0(path2analysis,
                                '/diff_peaks_7_3h_-IFN_polyC_VS_5_3h_pIFN_polyC.tsv.gz'))


## DEGs comparing 3 hrs dsRNA vs polyC during IFN conditions
if ( rerun_analysis == TRUE) {
        degs.stim <- FindMarkers(
                object = subset(seurat,
                                samples %in% c('5_3h_pIFN_polyC', '4_3h_pIFN_dsRNA')),
                ident.1 = '4_3h_pIFN_dsRNA',
                ident.2 = '5_3h_pIFN_polyC',
                min.pct = 0,
                logfc.threshold = 0,
                test.use = "wilcox"
        )
        colnames(degs.stim) <- paste0(colnames(degs.stim), '_pVsndsRNA_IFN')
        degs.stim <- rownames_to_column(degs.stim, var = 'peak')
        write_tsv(degs.stim,
                  file = paste0(path2analysis,
                                '/diff_peaks_4_3h_pIFN_dsRNA_VS_5_3h_pIFN_polyC.tsv.gz'))
}

degs.stim  <- read_tsv(paste0(path2analysis,
                              '/diff_peaks_4_3h_pIFN_dsRNA_VS_5_3h_pIFN_polyC.tsv.gz'))


## Massaging results to make a vulcano plot
degs.priming <- read_tsv(file = paste0(path2analysis,
                                       '/diff_peaks_7_3h_-IFN_polyC_VS_5_3h_pIFN_polyC.tsv.gz')) 
degs.stim <- read_tsv(file = paste0(path2analysis,
                                       '/diff_peaks_4_3h_pIFN_dsRNA_VS_5_3h_pIFN_polyC.tsv.gz'))
colnames(degs.priming); colnames(degs.stim)
peaks <- merge(degs.priming, degs.stim)
head(peaks)



## checking direction of peak regulation
top.peaks <- peaks %>%
        filter(p_val_pVsnIFN_polyC < 0.05) %>%
        arrange(desc(avg_log2FC_pVsnIFN_polyC)) %>%
        head() %>%
        pull(peak)

dotplot.check <- DotPlot(subset(seurat,
               samples %in% c('7_3h_-IFN_polyC', '5_3h_pIFN_polyC',
                              '4_3h_pIFN_dsRNA')), 
        group.by = 'samples',
        features = top.peaks,
        cols = 'RdYlGn') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        coord_flip() +
        labs(x='', y='')


head(peaks)
##write_tsv(peaks,
##          '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/chromatin_priming_comparative_analysis/diff_peaks_priming_Vs_stim.tsv.gz')


## Plotting
peaks <- mutate(peaks, 
                label = ifelse(peak %in% top.peaks,
                               peak, '')) 
peaks <- mutate(peaks, highlight = ifelse(peak %in% top.peaks,
                                          TRUE, FALSE))
primingVsStim.scatter_plot <- peaks %>%
        ggplot(aes(x = avg_log2FC_pVsnIFN_polyC, 
                   y = avg_log2FC_pVsndsRNA_IFN,
                   label = label)) +
        geom_point(alpha=0.5,
                   colour = 'steelblue') +
        geom_text_repel(max.overlaps = 100, force = 5) +
        geom_point(data = filter(peaks,  highlight == TRUE),
                   aes(x = avg_log2FC_pVsnIFN_polyC, 
                       y = avg_log2FC_pVsndsRNA_IFN,
                       colour = 'red')) +
        geom_hline(yintercept = 0,
                   linetype = 'dashed',
                   colour = 'red') +
        geom_vline(xintercept = 0,
                   linetype = 'dashed',
                   colour = 'red') +
        geom_smooth(method = 'lm') +
        theme_classic() +
        theme(legend.position = 'none') +
        labs(x='Log2FC IFN (+/-) | polyC',
             y= 'Log2FC dsRNA(+/-) | IFN')


primingVsStim.scatter_plot + dotplot.check
  

pdf(file = paste0(path2figures, '/sattter_plot_priming_vs_stim.pdf'),
    height = 6, width = 6)
primingVsStim.scatter_plot
dev.off()


pdf(file = paste0(path2figures, '/dotplot_accesibility_selected_markers_check.pdf'),
    height = 5, width = 6)
dotplot.check
dev.off()



##################################
## labeling peaks according to up- or down-regulation
peaks <- peaks %>%
        mutate(priming=ifelse(avg_log2FC_pVsnIFN_polyC>=0,
                              'primed', 'non-primed')) %>%
        mutate(stimulation=ifelse(avg_log2FC_pVsndsRNA_IFN>=0,
                                  'stimulated', 'non-stimulated'))



bar.perc <- peaks %>%
        ggplot(aes(x=stimulation,
                   fill=priming)) +
                geom_bar(position = 'fill') +
        scale_y_continuous(labels = scales::percent) +
        scale_fill_manual(values = priming_categories) +
        labs(x='', y='% Peaks',
             fill='') +
        coord_flip() +
        theme_bw() +
        theme(legend.position = 'bottom')

bar.numb <- peaks %>%
        ggplot(aes(x=stimulation,
                   fill=priming)) +
        geom_bar(position = 'stack') +
        scale_fill_manual(values = priming_categories) +
        labs(x='', y='# Peaks',
             fill='') +
        coord_flip() +
        theme_bw() +
        theme(legend.position = 'bottom')



pdf(file = paste0(path2figures, '/barplot_primingVsStimulated_peaks.pdf'),
    height = 5, width = 8)
bar.perc / bar.numb
dev.off()


##################################
## Combining Chromatin and Gene expression information
seurat.gex <- readRDS('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/seurat/ifn_treated_cells_rnaseq_processed_seurat.rds')

seurat.gex$treatment %>% 
        table()

## Calculating degs after stimulation
Idents(seurat.gex) <- seurat.gex$treatment 
if ( rerun_analysis == TRUE ) {
        markers.stim <- FindMarkers(
                subset( seurat.gex, treatment %in% c('3h_pIFN_polyC',
                                                     '3h_pIFN_dsRNA') ), 
                ident.1 = '3h_pIFN_polyC', 
                ident.2 = '3h_pIFN_dsRNA', 
                logfc.threshold = 0, 
                min.pct = 0
        )
        markers.stim <- rownames_to_column(markers.stim, var = 'gene_annotation')
        write_tsv(markers.stim,
                  paste0(path2analysis, 'degs_priming_3h_pIFN_dsRNA_VS_3h_pIFN_polyC.tsv.gz'))
}

markers.stim <- read_tsv(paste0(path2analysis, 'degs_priming_3h_pIFN_dsRNA_VS_3h_pIFN_polyC.tsv.gz'))

## Loading annotated peaks to genes
peaks.ann <- read_rds('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/peaks_annotation/annotated_peaks_unambiguously.rds')
peaks.ann <- mutate(peaks.ann, peak=paste0(seqnames, ':', start, '-', end))
head(peaks.ann)
length(unique(peaks.ann$peak))==length(peaks.ann$peak)
peaks$'gene_annotation' <- plyr::mapvalues(peaks$peak, 
                                           from = peaks.ann$peak,
                                           to = peaks.ann$annot.symbol) 


intersect(colnames(peaks), colnames(markers.stim))
peaks_gex <- merge(peaks, markers.stim)
head(peaks_gex)


## Subsetting to IFN genes
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



peaks_gex %>%
        filter(p_val < 0.05 & p_val_pVsnIFN_polyC < 0.05) %>%
        ggplot(aes(x=avg_log2FC_pVsnIFN_polyC,
                   y=avg_log2FC)) +
                geom_point()

peaks_gex %>%
        ggplot(aes(x = avg_log2FC_pVsnIFN_polyC, 
                   y = avg_log2FC_pVsndsRNA_IFN,
                   label = label)) +
        geom_point(alpha=0.5,
                   colour = 'steelblue') +
        #geom_text_repel(max.overlaps = 100, force = 5) +
        geom_point(data = filter(peaks_gex,
                                 gene_annotation %in% geneSets),
                   aes(x = avg_log2FC_pVsnIFN_polyC, 
                       y = avg_log2FC_pVsndsRNA_IFN,
                       colour = 'red'),
                   alpha = 0.5 ) +
        geom_hline(yintercept = 0,
                   linetype = 'dashed',
                   colour = 'red') +
        geom_vline(xintercept = 0,
                   linetype = 'dashed',
                   colour = 'red') +
        geom_smooth(method = 'lm') +
        theme_classic() +
        theme(legend.position = 'none') +
        labs(x='Log2FC IFN (+/-) | polyC',
             y= 'Log2FC dsRNA(+/-) | IFN')


