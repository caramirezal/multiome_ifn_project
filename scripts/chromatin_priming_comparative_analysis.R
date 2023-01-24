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
library(enrichR)

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


## Loading annotated peaks to genes
peaks.ann <- read_rds('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/analysis/peaks_annotation/annotated_peaks_unambiguously.rds')
peaks.ann <- mutate(peaks.ann, peak=paste0(seqnames, ':', start, '-', end))
head(peaks.ann)
length(unique(peaks.ann$peak))==length(peaks.ann$peak)




###############################
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





#######################################################
## Differential peaks comparing 3 hrs IFN vs non-IFN (during polyC)
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





#############################################################
## Differential peaks comparing 3 hrs dsRNA vs polyC during IFN conditions
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

degs.priming %>%
        mutate(gene_annotation=plyr::mapvalues(degs.priming$peak,
                                               from = peaks.ann$peak,
                                               to = peaks.ann$annot.symbol),
               annot.type=plyr::mapvalues(degs.priming$peak,
                                               from = peaks.ann$peak,
                                               to = peaks.ann$annot.type)) %>%
        filter(p_val_pVsnIFN_polyC < 0.05 & annot.type == 'hg38_genes_promoters') %>%
        arrange(desc(avg_log2FC_pVsnIFN_polyC)) %>%
        as.data.frame() %>%
        head()
degs.stim <- read_tsv(file = paste0(path2analysis,
                                       '/diff_peaks_4_3h_pIFN_dsRNA_VS_5_3h_pIFN_polyC.tsv.gz'))
degs.stim %>%
        mutate(gene_annotation=plyr::mapvalues(degs.stim$peak,
                                               from = peaks.ann$peak,
                                               to = peaks.ann$annot.symbol),
               annot.type=plyr::mapvalues(degs.stim$peak,
                                          from = peaks.ann$peak,
                                          to = peaks.ann$annot.type)) %>%
        filter(p_val_pVsndsRNA_IFN < 0.05 & annot.type == 'hg38_genes_promoters') %>%
        arrange(desc(avg_log2FC_pVsndsRNA_IFN)) %>%
        as.data.frame() %>%
        head()
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




####################################################
## Plotting priming Vs stimulation
peaks <- mutate(peaks, 
                label = ifelse(peak %in% top.peaks,
                               peak, '')) 
peaks <- mutate(peaks, highlight = ifelse(peak %in% top.peaks,
                                          TRUE, FALSE))
primingVsStim.scatter_plot <- peaks %>%
        filter(`pct.1_pVsnIFN_polyC` > 0.05 & 
                       `pct.1_pVsndsRNA_IFN` > 0.05 ) %>%
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
        theme(legend.position = 'none',
              text = element_text(size = 20)) +
        labs(x='Log2FC IFN (+/-) | polyC',
             y= 'Log2FC dsRNA(+/-) | IFN')




##################################################
## Plotting frequency of peaks with a certain % of cells with counts > 0
 pct.cells.hist <- select(peaks, `pct.1_pVsnIFN_polyC`, `pct.1_pVsndsRNA_IFN`) %>%
        reshape2::melt() %>%
        ggplot(aes(x=100*value)) +
                geom_histogram(bins = 50) +
                geom_vline(xintercept = 5,
                           linetype = 'dashed',
                           colour = 'red') +
                facet_wrap(~variable, ncol=1) +
                scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
                theme_bw() +
                theme(text = element_text(size=20)) +
                labs(x='% Cells with counts > 0', y='')


primingVsStim.scatter_plot + dotplot.check
  

pdf(file = paste0(path2figures, '/scatter_plot_priming_vs_stim.pdf'),
    height = 6, width = 6)
primingVsStim.scatter_plot
pct.cells.hist
dev.off()


pdf(file = paste0(path2figures, '/dotplot_accesibility_selected_markers_check.pdf'),
    height = 5, width = 6)
dotplot.check
dev.off()





#################################################
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




##################################################
## Combining Chromatin and Gene expression information
seurat.gex <- readRDS('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/seurat/ifn_treated_cells_rnaseq_processed_seurat.rds')

seurat.gex$treatment %>% 
        table()



#################################################
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
                  paste0(path2analysis, '/degs_priming_3h_pIFN_dsRNA_VS_3h_pIFN_polyC.tsv.gz'))
}

markers.stim <- read_tsv(paste0(path2analysis, '/degs_priming_3h_pIFN_dsRNA_VS_3h_pIFN_polyC.tsv.gz'))

        


######################################
## Comparing transcriptional signatures

## Calculating degs after stimulation
Idents(seurat.gex) <- seurat.gex$treatment 
if ( rerun_analysis == TRUE ) {
        markers.prim <- FindMarkers(
                subset( seurat.gex, treatment %in% c('3h_pIFN_polyC',
                                                     '3h_-IFN_polyC') ), 
                ident.1 = '3h_pIFN_polyC', 
                ident.2 = '3h_-IFN_polyC', 
                logfc.threshold = 0, 
                min.pct = 0
        )
        markers.prim <- rownames_to_column(markers.prim, var = 'gene_annotation')
        write_tsv(markers.prim,
                  paste0(path2analysis, '/degs_priming_3h_pIFN_polyCA_VS_3h_-IFN_polyC.tsv.gz'))
}

markers.prim <- read_tsv(paste0(path2analysis, '/degs_priming_3h_pIFN_polyCA_VS_3h_-IFN_polyC.tsv.gz'))


colnames(markers.prim) <- c('gene', paste0(colnames(markers.prim), '_priming')[2:ncol(markers.prim)] )
colnames(markers.stim) <- c('gene', paste0(colnames(markers.stim), '_stim')[2:ncol(markers.stim)])
markers <- merge(markers.prim, markers.stim)
head(markers)
markers <- markers %>%
        mutate(highlight=ifelse(gene %in% geneSets,
                                TRUE, FALSE))



#####################################
## GEX Priming signature
up.reg.prim <- arrange(markers.prim, desc(avg_log2FC_priming)) %>%
                        head(5) %>%
                        pull(gene)
down.reg.prim <- arrange(markers.prim, desc(avg_log2FC_priming)) %>%
        tail(5) %>%
        pull(gene)
markers.prim <- markers.prim %>%
        mutate(top=ifelse(gene %in% c(up.reg.prim,
                                      down.reg.prim),
                          TRUE, FALSE) ) %>%
        mutate(label=ifelse(top==TRUE, gene, ''))
vulcano.prim.gex <- markers.prim %>%
        ggplot(aes(x=avg_log2FC_priming,
                               y=-log10(p_val_priming),
                               label=label)) +
        geom_point(colour='gray80') +
        geom_point(data = subset(markers.prim,
                                 top == TRUE),
                   aes(x=avg_log2FC_priming,
                       y=-log10(p_val_priming)),
                   colour = 'red') +
        geom_text_repel(max.overlaps = 100) +
        geom_vline(xintercept = 0, 
                   linetype = 'dashed',
                   colour = 'red') +
        theme_bw() +
        theme(text = element_text(size=20)) +
        labs(x='Avg Log2FC',
             y='-Log10(p-Value)',
             title = 'IFN (+/-) | polyC') +
        xlim(c(-7, 7))



## check primed genes
dotplot.prim <- DotPlot(object = subset(seurat.gex, 
                        treatment %in% c('3h_pIFN_polyC',
                                         '3h_-IFN_polyC')), 
        features = c(up.reg.prim, down.reg.prim)) +
        scale_colour_viridis() +
        coord_flip() +
        labs(x='', y='') +
                theme(axis.text.x = element_text(hjust = 1,
                                                 angle=45),
                      text = element_text(size=20))


dotplot.stim <- DotPlot(object = subset(seurat.gex, 
                                        treatment %in% c('3h_pIFN_polyC',
                                                         '3h_pIFN_dsRNA')), 
                        features = c(up.reg.stim, down.reg.stim)) +
        scale_colour_viridis() +
        coord_flip() +
        labs(x='', y='') +
        theme(axis.text.x = element_text(hjust = 1,
                                         angle=45),
              text = element_text(size=20))


pdf(paste0(path2figures, '/dotplot_checks_gex_signatures.pdf'),
    width = 10, height = 5)
dotplot.prim + dotplot.stim
dev.off()

#####################################
## GEX Stimulation signature
up.reg.stim <- arrange(markers.stim, desc(avg_log2FC_stim)) %>%
        head(5) %>%
        pull(gene)
down.reg.stim <- arrange(markers.stim, desc(avg_log2FC_stim)) %>%
        tail(5) %>%
        pull(gene)
markers.stim <- markers.stim %>%
        mutate(top=ifelse(gene %in% c(up.reg.stim,
                                      down.reg.stim),
                          TRUE, FALSE) ) %>%
        mutate(label=ifelse(top==TRUE, gene, ''))
markers.stim <- mutate(markers.stim, avg_log2FC_stim=-avg_log2FC_stim)      ## 
vulcano.stim.gex <- markers.stim %>%
        ggplot(aes(x=avg_log2FC_stim,
                   y=-log10(p_val_stim),
                   label=label)) +
        geom_point(colour='gray80') +
        geom_point(data = subset(markers.stim,
                                 top == TRUE),
                   aes(x=avg_log2FC_stim,
                       y=-log10(p_val_stim)),
                   colour = 'red') +
        geom_text_repel(max.overlaps = 100) +
        geom_vline(xintercept = 0, 
                   linetype = 'dashed',
                   colour = 'red') +
        theme_bw() +
        theme(text = element_text(size=20)) +
        labs(x='Avg Log2FC',
             y='-Log10(p-Value)',
             title = 'dsRNA (+/-) | IFN+') +
        xlim(c(-7, 7))



pdf(paste0(path2figures, '/vulcano_plots_priming_stim_signatures_gex.pdf'),
    width = 6, height = 9)
vulcano.prim.gex / vulcano.stim.gex
dev.off()


        
        
pdf(paste0(path2figures, '/scatter_plot_priming_Vs_stim_signs_gex.pdf'),
    width = 5, height = 5)
markers %>%        
        filter( pct.1_stim > 0.1 & pct.1_priming > 0.1 ) %>%
        ggplot(aes(x=avg_log2FC_priming,
                   y=-avg_log2FC_stim)) +
                geom_point(colour='gray80') +
                geom_point(data = subset(markers, 
                                         highlight == TRUE),
                           aes(x=avg_log2FC_priming,
                               y=avg_log2FC_stim),
                           colour='red') +
                geom_hline(yintercept = 0,
                           linetype = 'dashed',
                           colour = 'black') +
                geom_vline(xintercept = 0,
                           linetype = 'dashed',
                           colour = 'black') +
                geom_smooth(method = 'lm') +
                scale_color_viridis() +
                theme_bw() +
                theme(text = element_text(size=20)) +
                labs(x='Log2FC IFN (+/-) | polyC',
                     y='Log2FC dsRNA (+/-) | IFN') 
dev.off()




#################################################
## Calculating genes after priming




peaks$'gene_annotation' <- plyr::mapvalues(peaks$peak, 
                                           from = peaks.ann$peak,
                                           to = peaks.ann$annot.symbol) 
peaks$'annot.type' <- plyr::mapvalues(peaks$peak, 
                                      from = peaks.ann$peak,
                                      to = peaks.ann$annot.type) 

intersect(colnames(peaks), colnames(markers.stim))
peaks_gex <- merge(peaks, markers.stim)
head(peaks_gex)




peaks_gex %>%
        write_tsv(file = paste0(path2analysis, '/diff_peaks_gex_ann.tsv.gz'))




######################################
## Split the genes in three categories:
## a) gained after IFN
## b) gained after dsRNA
peaks_gex <- peaks_gex %>%
        mutate(priming=ifelse(avg_log2FC_pVsnIFN_polyC>=0 & p_val_pVsnIFN_polyC < 0.05,
                              'gained_after_IFN', 'non-gained')) %>%
        mutate(stimulation=ifelse(avg_log2FC_pVsndsRNA_IFN>=0 & p_val_pVsndsRNA_IFN < 0.05,
                                  'gained_after_dsRNA', 'non-gained')) %>%
        mutate(gex_upregulation=ifelse(avg_log2FC>=0 & p_val < 0.05 ,
                                       'up-regulated', 'non-upregulated')) %>%
        mutate(category=paste0(priming, ':', stimulation))
peaks_gex$category <- plyr::mapvalues(peaks_gex$category, 
                                        from = unique(peaks_gex$category),
                                        to = c('non-gained',
                                               'gained_after_IFN',
                                               'gained_after_dsRNA'))






######################################
## Priming Signature - Vulcano plot
primed.peaks.df <- subset(peaks_gex, 
                          avg_log2FC_pVsnIFN_polyC > 0 & 
                                  p_val_pVsnIFN_polyC < 0.05 &
                                  annot.type == 'hg38_genes_promoters') %>%
                        arrange(desc(avg_log2FC_pVsnIFN_polyC))
nb.priming_peaks <- primed.peaks.df %>% nrow()
selected.priming.peaks <- head(primed.peaks.df) 
top3.downregulated <- peaks_gex %>%
                           filter(p_val_pVsnIFN_polyC < 0.05 ) %>%
                           arrange(desc(avg_log2FC_pVsnIFN_polyC)) %>%
                                   tail(n=3) %>%
                           pull(peak) 
vulcano.priming <- peaks_gex %>%
        mutate(label=ifelse(peak %in% c(primed.peaks.df$peak[1:3],
                                        top3.downregulated),  
                            gene_annotation, '')) %>%
        ggplot(aes(x=avg_log2FC_pVsnIFN_polyC,
                   y=-log10(p_val_pVsnIFN_polyC),
                   label=label)) +
        geom_point(colour='gray80') +
        geom_point(data = subset(peaks_gex, 
                                 avg_log2FC_pVsnIFN_polyC > 0 & 
                                         p_val_pVsnIFN_polyC < 0.05 &
                                         annot.type == 'hg38_genes_promoters'),
                   aes(x=avg_log2FC_pVsnIFN_polyC,
                       y=-log10(p_val_pVsnIFN_polyC)),
                   colour = 'darkorange',
                   alpha=0.5) +
        geom_text(data = data.frame(
                'avg_log2FC_pVsnIFN_polyC' = 0.1,
                'p_val_pVsnIFN_polyC' = 0.001
        ), aes(x=avg_log2FC_pVsnIFN_polyC,
               y=-log10(p_val_pVsnIFN_polyC),
               label = nb.priming_peaks),
        size = 12
        ) +
        geom_text_repel() +
        theme_classic() + 
        theme(axis.text = element_text(size=20),
              text = element_text(size=24)) +
        geom_vline(xintercept = 0,
                   linetype = 'dashed',
                   colour = 'red')  + 
        geom_hline(yintercept = -log10(0.05),
                   linetype = 'dashed',
                   colour = 'red') +
        xlim(c(-0.2, 0.2)) +
        labs(x='Avg Log2FC',
             y='-Log10(p-Value)',
             title = 'IFN (+/-) | polyC')






######################################
## Stimulation Signature - Vulcano plot
stim.peaks.df <- subset(peaks_gex, 
                          avg_log2FC_pVsndsRNA_IFN > 0 & 
                                  p_val_pVsndsRNA_IFN < 0.05 &
                                  annot.type == 'hg38_genes_promoters') %>%
        arrange(desc(avg_log2FC_pVsndsRNA_IFN))
nb.stim_peaks <- stim.peaks.df %>% nrow()
head(stim.peaks.df) 
top3.downregulated <- peaks_gex %>%
        filter(p_val_pVsndsRNA_IFN < 0.05 ) %>%
        arrange(desc(avg_log2FC_pVsndsRNA_IFN)) %>%
        tail(n=3) %>%
        pull(peak) 
vulcano.stim <- peaks_gex %>%
        mutate(label=ifelse(peak %in% c(stim.peaks.df$peak[1:3],
                                        top3.downregulated),  
                            gene_annotation, '')) %>%
        ggplot(aes(x=avg_log2FC_pVsndsRNA_IFN,
                   y=-log10(p_val_pVsndsRNA_IFN),
                   label=label)) +
        geom_point(colour='gray80') +
        geom_point(data = subset(peaks_gex, 
                                         avg_log2FC_pVsndsRNA_IFN > 0 & 
                                         p_val_pVsndsRNA_IFN < 0.05 &
                                         annot.type == 'hg38_genes_promoters'),
                   aes(x=avg_log2FC_pVsndsRNA_IFN,
                       y=-log10(p_val_pVsndsRNA_IFN)),
                   colour = 'darkolivegreen3',
                   alpha=0.5) +
        geom_text(data = data.frame(
                'avg_log2FC_pVsndsRNA_IFN' = 0.1,
                'p_val_pVsndsRNA_IFN' = 0.001
        ), aes(x=avg_log2FC_pVsndsRNA_IFN,
               y=-log10(p_val_pVsndsRNA_IFN),
               label = nb.stim_peaks),
        size = 12
        ) +
        geom_text_repel() +
        theme_classic() + 
        theme(axis.text = element_text(size=20),
              text = element_text(size=24)) +
        geom_vline(xintercept = 0,
                   linetype = 'dashed',
                   colour = 'red')  + 
        geom_hline(yintercept = -log10(0.05),
                   linetype = 'dashed',
                   colour = 'red') +
        xlim(c(-0.2, 0.2)) +
        labs(x='Avg Log2FC',
             y='-Log10(p-Value)',
             title = 'dsRNA (+/-) | IFN+')




###########################################
## Markers dsRNA stimulation
overlap.df <- rbind(
        subset(primed.peaks.df, p_val < 0.05 ),
        subset(stim.peaks.df, p_val < 0.05)
) %>% arrange(desc(avg_log2FC)) 
vulcano.stim.degs <- peaks_gex %>%
        ggplot(aes(x=avg_log2FC,
                   y=-log10(p_val))) +
        geom_point(colour='gray80') +
        ## Add primed genes
        geom_point(data = subset(peaks_gex, 
                                 p_val < 0.05 & 
                                 avg_log2FC_pVsndsRNA_IFN > 0 & 
                                         p_val_pVsndsRNA_IFN < 0.05 &
                                 annot.type == 'hg38_genes_promoters'),
                   aes(x=avg_log2FC,
                       y=-log10(p_val)),
                   colour = 'darkolivegreen3',
                   alpha=0.5) +
        ## Add stimulated genes
        geom_point(data = subset(peaks_gex, 
                                 p_val < 0.05 & 
                                         avg_log2FC_pVsnIFN_polyC > 0 & 
                                         p_val_pVsnIFN_polyC < 0.05 &
                                         annot.type == 'hg38_genes_promoters'),
                   aes(x=avg_log2FC,
                       y=-log10(p_val)),
                   colour = 'darkorange',
                   alpha=0.5) +
        ## up-regulated genes
        geom_text(data = data.frame(
                'avg_log2FC' = 2,
                'p_val' = (1/10**100)
        ), aes(x=avg_log2FC,
               y=-log10(p_val),
               label = nrow(subset(overlap.df, 0 < avg_log2FC))),
        size = 12
        ) +
        ## down-regulated genes
        geom_text(data = data.frame(
                'avg_log2FC' = -2.5,
                'p_val' = (1/10**100)
        ), aes(x=avg_log2FC,
               y=-log10(p_val),
               label = nrow(subset(overlap.df, avg_log2FC < 0 ))),
        size = 12
        ) +
        ## Add synergistic genes
        geom_vline(xintercept = 0,
                   linetype = 'dashed',
                   colour = 'red') +
        theme_classic() +
        theme(axis.text = element_text(size=20),
              text = element_text(size=24)) +
        labs(x='Avg Log2FC',
             y='-Log10(p-Value)',
             title = 'DEG - dsRNA (+/-) | IFN+')
                        


primed.genes <- subset(primed.peaks.df, p_val < 0.05 & 0 <  avg_log2FC ) %>% pull(gene_annotation)
stim.genes <- subset( stim.peaks.df, p_val < 0.05 & 0 <  avg_log2FC  ) %>% pull(gene_annotation)
intersect.genes <- intersect(primed.genes, stim.genes)
overlap.df$'condition_of_increase' <- sapply(
        1:nrow(overlap.df),
        function(idx){
              res <- 'None'
              if ( overlap.df$gene_annotation[idx] %in% stim.genes ) {
                      res <- 'dsRNA'
              } 
              if ( overlap.df$gene_annotation[idx] %in% primed.genes) {
                      res <- 'IFN'
              }
              if ( overlap.df$gene_annotation[idx] %in% stim.genes &
                   overlap.df$gene_annotation[idx] %in% primed.genes ) {
                      res <- 'Both'
              }
              res
        }
) 
barplot.primVsStim <- overlap.df %>%
        filter(avg_log2FC > 0 ) %>%
        ggplot(aes(x=gex_upregulation,
                   fill=condition_of_increase)) +
        geom_bar(position = 'stack') +
        scale_fill_manual(values = c('IFN'='darkorange',
                                     'dsRNA'='darkolivegreen3',
                                     'Both'='cyan3')) +
        theme_classic() +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.text = element_text(size=20),
              text= element_text(size=24),
              legend.position = 'bottom') +
        labs(x='', 
             y='# Up-regulated genes',
             fill='') +
        coord_flip()


pdf(paste0(path2figures, '/vulcano_priming_and_stim.pdf'),
    height = 8, width = 5)
( vulcano.priming / vulcano.stim ) 
dev.off()


pdf(paste0(path2figures, '/vulcano_degs_stimulation.pdf'),
    height = 7, width = 7)
vulcano.stim.degs 
dev.off()


pdf(paste0(path2figures, '/barplot_upregulated_degs_with_associated_peaks.pdf'),
    height = 2, width = 6)
barplot.primVsStim 
dev.off()




#############################################
## Gene set enrich analysis of genes up-regulated and with changes
## in chromatin accessibility
dbs <- listEnrichrDbs()
dbs <- c('Reactome_2022', 'GO_Biological_Process_2021', 'GO_Cellular_Component_2021',
         'KEGG_2021_Human', 'MSigDB_Hallmark_2020')
up_regulated_genes <- filter(primed.peaks.df, avg_log2FC > 0 &
                                     p_val < 0.05 ) %>% pull(gene_annotation)

enriched <- enrichr(up_regulated_genes, dbs)
enrich.plot <- plotEnrich(enriched$Reactome_2022, 
           showTerms = 20, 
           numChar = 60, 
           y = "Count", 
           orderBy = "P.value") +
                labs(x='', 
                     title='Up-regulated After IFN',
                     subtitle = 'Open after ')


enriched.primed <- enrichr(primed.genes, dbs)
enrich.primed.plot <- plotEnrich(enriched.primed$Reactome_2022, 
                          showTerms = 20, 
                          numChar = 60, 
                          y = "Count", 
                          orderBy = "P.value") +
        labs(x='', 
             title='Up-regulated After IFN')

enriched.stim <- enrichr(stim.genes, dbs)
enrich.stim.plot <- plotEnrich(enriched.stim$Reactome_2022, 
                          showTerms = 20, 
                          numChar = 60, 
                          y = "Count", 
                          orderBy = "P.value") +
        labs(x='', 
             title='Up-regulated After dsRNA')


pdf(paste0(path2figures, '/gsea_degs_primed_and_stimulated.pdf'),
    width = 20, height = 5)
enrich.primed.plot + enrich.stim.plot
dev.off()

intersect(primed.genes, geneSets)





################################################
## Checking the state of primed open regions after dsRNA stimulation
write_tsv(primed.peaks.df,
          file = paste0(path2analysis, '/genes_upregualted_after_priming_487.tsv.gz'))
write_tsv(stim.peaks.df, 
          file = paste0(path2analysis, '/genes_upregulated_after_stimulation_219.tsv.gz'))


primed.peaks.df <- read_tsv(paste0(path2analysis, '/genes_upregualted_after_priming_487.tsv.gz'))
primed.peaks.df <- mutate(primed.peaks.df,
                          category = sapply(1:nrow(primed.peaks.df),
                                            function(idx){
                                                    res <- 'non_significant_change'
                                                    if ( primed.peaks.df$avg_log2FC_pVsndsRNA_IFN[idx] > 0 &
                                                         primed.peaks.df$p_val_pVsndsRNA_IFN[idx] < 0.05 ) {
                                                            res <- 'upregulated_after_stimulation'
                                                    } 
                                                    if ( primed.peaks.df$avg_log2FC_pVsndsRNA_IFN[idx] < 0 &
                                                         primed.peaks.df$p_val_pVsndsRNA_IFN[idx] < 0.05 ) {
                                                            res <- 'down_after_stimulation'
                                                    }
                                                    return(res)
                                            }))
primed.peaks.df$category %>% table()


#######################################
## Visualization of the peaks
pdf(paste0(path2figures, '/projected_primed_peaks_on_stimulation_vulcano.pdf'),
    height = 4, width = )
vulcano.priming +
peaks_gex %>%
        ggplot(aes(x=avg_log2FC_pVsndsRNA_IFN,
                   y=-log10(p_val_pVsndsRNA_IFN))) +
        geom_point(colour='gray80') +
        geom_point(data = subset(peaks_gex, 
                                 peak %in% primed.peaks.df$peak),
                   aes(x=avg_log2FC_pVsndsRNA_IFN,
                       y=-log10(p_val_pVsndsRNA_IFN)),
                   colour = 'darkorange',
                   alpha=0.5) +
        theme_classic() + 
        theme(axis.text = element_text(size=20),
              text = element_text(size=24)) +
        geom_vline(xintercept = 0,
                   linetype = 'dashed',
                   colour = 'red')  + 
        geom_hline(yintercept = -log10(0.05),
                   linetype = 'dashed',
                   colour = 'red') +
        xlim(c(-0.2, 0.2)) +
        labs(x='Avg Log2FC',
             y='-Log10(p-Value)',
             title = 'dsRNA (+/-) | IFN+')
dev.off()





