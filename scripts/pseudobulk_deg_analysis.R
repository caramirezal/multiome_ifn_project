## This script contains a DE analysis of human cell lines pre-treated (or not) with IFN
## and stimulated (or not) with dsRNA

library(Seurat)
library(edgeR)
library(scater)
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(viridis)


## loading settings
source('/media/ag-cherrmann/cramirez/multiome_ifn_project/scripts/settings.R')

## Loading data in seurat format
seurat <- readRDS('/media/ag-cherrmann/cramirez/multiome_ifn_project/data/seurat/ifn_treated_cells_rnaseq_processed_seurat.rds')
## checking data
#DimPlot(seurat, group.by = 'orig.ident', cols = sample.colors)

## creating folder to store results
path2figures <- '/media/ag-cherrmann/cramirez/multiome_ifn_project/figures/pseudobulk_deg_analysis'
if ( ! dir.exists(path2figures)){
    dir.create(path2figures)
}

## creating folder to store analysis
path2analysis <- '/media/ag-cherrmann/cramirez/multiome_ifn_project/analysis/pseudobulk_deg_analysis'
if ( ! dir.exists(path2analysis)){
    dir.create(path2analysis)
}

## Transforming seurat object to a single cell experiment object
scexp <- as.SingleCellExperiment(seurat)
colData(scexp)$'dsRNA' <- ifelse(grepl('dsRNA', seurat$stimulus), 'dsRNA', 'polyC') %>% as.factor()
colData(scexp)$'dsRNA' <- factor(colData(scexp)$'dsRNA', levels = c('polyC', 'dsRNA'))
colData(scexp)$'IFN' <- ifelse(grepl('pIFN', seurat$stimulus), '-IFN', 'pIFN') %>% as.factor()
colData(scexp)$'IFN' <- factor(colData(scexp)$'IFN', levels = c('-IFN', 'pIFN'))


## aggregating cells in samples samples
summed <- aggregateAcrossCells(scexp, 
                id = colData(scexp)[, c('ident')]
)

## creating dge list object
dge_list <- DGEList(counts(summed), samples=colData(summed))
dge_list

## filtering genes
keep <- filterByExpr(dge_list, group=summed$'ident')
dge_filtered <- dge_list[keep, ]


## Calculation of normalization factors
dge_filtered <- calcNormFactors(dge_filtered)


## Visualizing MD plots 
pdf(paste0(path2figures, '/mdplot_samples.pdf'),
    width = 12, height = 8
)
par(mfrow=c(2,3))
for ( i in seq_len(ncol(dge_filtered))){
    plotMD(dge_filtered, column=i,
            main=colnames(dge_filtered)[i]
            )   
}
dev.off()


## PCA visualization
pdf('/media/ag-cherrmann/cramirez/multiome_ifn_project/figures/pseudobulk_deg_analysis/pca_mds.pdf')
limma::plotMDS(edgeR::cpm(dge_filtered, log=TRUE), cex=0.8)
dev.off()


## creating the design matrix
design <- model.matrix(~factor(dsRNA) + factor(IFN, levels = c('pIFN', '-IFN')) + factor(dsRNA)*factor(IFN, levels = c('pIFN', '-IFN')), 
                                dge_filtered$samples)


## Evaluation of the mean to variance distribution
dge_filtered <- estimateDisp(dge_filtered, design)


pdf(paste0(path2figures, '/BCV_plot.pdf'))
plotBCV(dge_filtered)
dev.off()


## Using a Quasi-likehood for mean-variance dispersion
design <- model.matrix(~factor(dsRNA) + factor(IFN, levels = c('pIFN', '-IFN')) + factor(dsRNA)*factor(IFN, levels = c('pIFN', '-IFN')), 
                                dge_filtered$samples)
fit <- glmQLFit(dge_filtered, design, robust=TRUE)
fit$coefficients %>%
    head()
fit.df <- as.data.frame(fit$coefficient)
fit.df$gene <- rownames(fit.df)
write_tsv(fit.df,
    file = '/media/ag-cherrmann/cramirez/multiome_ifn_project/analysis/pseudobulk_deg_analysis/coefficients_edgr_model.tsv.gz'
)


## performing a likehood ratio test
lrt <- glmLRT(fit, method='fdr')
pvals.df <- lrt$table %>% as.data.frame()
pvals.df$'gene' <- rownames(pvals.df)
head(pvals.df)
fit.df <- merge(fit.df, pvals.df)
## filtering p-vals > 0.05
fit.df <- filter(fit.df, PValue< 0.05)

top5.dsRNA <- fit.df %>%
            arrange(desc(`factor(dsRNA)dsRNA`)) %>%
                head(5) %>% pull(gene) 
top5.ifn <- fit.df %>%
            arrange(desc(`factor(IFN, levels = c("pIFN", "-IFN"))-IFN`)) %>%
                head(5) %>% pull(gene)
top5.interact <- fit.df %>%
            arrange(desc(`factor(dsRNA)dsRNA:factor(IFN, levels = c("pIFN", "-IFN"))-IFN`)) %>%
                head(5) %>% pull(gene)
tail5.dsRNA <- fit.df %>%
            arrange(desc(`factor(dsRNA)dsRNA`)) %>%
                tail(5) %>% pull(gene)
tail5.ifn <- fit.df %>%
            arrange(desc(`factor(IFN, levels = c("pIFN", "-IFN"))-IFN`)) %>%
                tail(5) %>% pull(gene)

pdf('/media/ag-cherrmann/cramirez/multiome_ifn_project/figures/pseudobulk_deg_analysis/coefficients_edgr_model.pdf',
    width=8, height=8)
fit.df %>% 
    filter(PValue< 0.05) %>%
    mutate(label=ifelse(gene %in% c(top5.dsRNA, top5.ifn, top5.interact,
                                    tail5.dsRNA, tail5.ifn),
                        gene, '')
                        )  %>%
    ggplot(aes(x=`factor(dsRNA)dsRNA`, 
               y=`factor(IFN, levels = c("pIFN", "-IFN"))-IFN`,
               colour=`factor(dsRNA)dsRNA:factor(IFN, levels = c("pIFN", "-IFN"))-IFN`,
               label=label)
) +
    geom_point() +
    geom_hline(yintercept=0, linetype='dashed',
                colour='red') +
    geom_vline(xintercept=0, linetype='dashed',
                colour='red') +
    geom_smooth(method='lm') +
    stat_regline_equation(label.y = 6, 
                          aes(label = ..eq.label..)) +
    stat_regline_equation(label.y = 5, 
                          aes(label = ..rr.label..)) +
    geom_text_repel(max.overlaps=10000,
                    force=5) +
    xlim(-10, 10) +
    ylim(-10, 10) +
    scale_colour_viridis() +
    theme_classic() +
    theme(panel.grid=element_blank()) +
    labs(x='dsRNA coefficient', y='IFN coefficient', colour='dsRNA*IFN coeff')
dev.off()


pdf(paste0(path2figures, '/QL_dispersion.pdf'))
plotQLDisp(fit)
dev.off()




## Inspection of top genes
summed <- logNormCounts(summed, size.factors=NULL)
pdf(paste0(path2figures, '/violinplot_deg_check_top5_ifn.pdf'),
    width=7, height=8)
plotExpression(summed,
                features=c(top5.ifn, top5.dsRNA, top5.interact),
                x='ident',
                ncol=5,
                colour_by=I(factor(summed$ident))) +
                theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank()
                      ) +
                      labs(x='') 

dev.off()