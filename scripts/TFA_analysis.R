## This script contain a transcription factor analysis of cells primed with IFN
## and stimulated with dsRNA
library(Seurat)
library(dplyr)
library(readr)
library(AUCell)
library(ComplexHeatmap)
library(ggrepel)
library(igraph)
library(viridis)

## initial settings
source('/media/ag-cherrmann/cramirez/multiome_ifn_project/scripts/settings.R')

## path to figures 
path2figures <- '/media/ag-cherrmann/cramirez/multiome_ifn_project/figures/tfa_analysis'
if ( ! dir.exists(path2figures) ) {
       dir.create(path2figures) 
}

## path to regulatory networks inferred with celloracle
path2networks <- '/media/ag-cherrmann/cramirez/multiome_ifn_project/analysis/celloracle/'
files <- list.files(path2networks, pattern = 'raw_GRN', full.names = TRUE)
names(files) <- gsub('.*raw_GRN_for_|\\.csv', '', files)
files


## getting transcription factor lists
all.tfs <- readLines('https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_curated_tfs.txt')
head(all.tfs)

##-------------------------------
networks.list <- lapply(files, read.csv)
lapply(networks.list, dim)
lapply(networks.list, function(x) x %>% pull(source) %>% unique %>% length)
lapply(networks.list, function(x) x %>% pull(target) %>% unique %>% length)
lapply(networks.list, head)


##-------------------------------
## Visualisation of the significances
pvals.list <- lapply(networks.list, 
                     function(df) {
                             p <- ggplot(df, aes(x=p)) +
                                        geom_histogram(bins = 100) +
                                     geom_vline(xintercept = 0.001, 
                                                linetype='dashed', 
                                                colour='red')
                             return(p)
                     })
gridExtra::grid.arrange(grobs=pvals.list, ncol=3)


##--------------------------------------------------------
## filtering networks by significances
networks.list <- lapply(networks.list, 
                   function(df){
                           dplyr::filter(df,
                                         p < 0.001) #%>%
                                   #arrange(desc(coef_mean)) %>%
                                   #head(5000)
                   })
lapply(networks.list, dim)
lapply(networks.list, function(x) x %>% pull(source) %>% unique %>% length)
lapply(networks.list, function(x) x %>% pull(target) %>% unique %>% length)
lapply(networks.list, head)



##----------------------------------------------------
## extracting regulons
conditions <- names(networks.list)
names(conditions) <- conditions
extract_regulon <- function(df){
        tfs <- df$source %>% unique() 
        names(tfs) <- tfs
        ## filtering only TFs
        tfs <- tfs[tfs %in% all.tfs]
        lapply(tfs, 
               function(tf) {
                       filter(df, source == tf) %>%
                               pull(target) 
               })
}
regulons.list.cond <- lapply(conditions,
                             function(cond){
                                     extract_regulon(networks.list[cond][[1]])
                             })



##--------------------------------
## Visualisation of the size of the regulons
tfs <- filter(networks.list$`3_16h_-IFN_dsRNA`,
              source %in% all.tfs) %>% pull(source) %>% unique()
names(tfs) <- tfs
regulon.sizes <- sapply(tfs,
       function(tf){
             filter(networks.list$`3_16h_-IFN_dsRNA`, source ==  tf) %>%
                       pull(target) %>%
                       unique() %>%
                       length()
       })
reg.sizes.df <- data.frame(tf=names(regulon.sizes),
                           regulon_size=regulon.sizes) 
reg.sizes.df %>%
        ggplot(aes(x=regulon_size)) +
                geom_histogram(bins = 100,
                               fill = 'steelblue') +
                geom_vline(xintercept = 30, linetype='dashed', colour='red') +
                theme_classic()
dim(reg.sizes.df)



## Extracting non-filtered network
regulons.non.filtered <- lapply(tfs,
                        function(tf){
                                filter(networks.list$`3_16h_-IFN_dsRNA`, 
                                       source ==  tf) %>%
                                        pull(target) %>%
                                        unique() 
                        })


##----------------------------------------------------
## filtering out regulons with less than 30 targets
keep.tf <- sapply(regulons.non.filtered, 
                                function(regulon) {
                                        length(regulon) >= 30 } )
reg.filt.bySize <- regulons.non.filtered[keep.tf]



##--------------------------------------------------
## Loading gene expression data
seurat <- readRDS('/media/ag-cherrmann/cramirez/multiome_ifn_project/data/seurat/ifn_treated_cells_rnaseq_processed_seurat.rds')



##--------------------------------------------------
## Calculation of Transcription factor activities
lapply(reg.filt.bySize, head)
## Calculating the scores
cells_rankings <- AUCell_buildRankings(seurat@assays$SCT@scale.data)
cells_AUC <- AUCell_calcAUC(reg.filt.bySize, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.05)
auc.mtx <- getAUC(cells_AUC)
auc.df <- as.data.frame(t(auc.mtx))
head(auc.df)



## Saving analysis
path2tfa_analysis <- '/media/ag-cherrmann/cramirez/multiome_ifn_project/analysis/tfa'
if ( ! dir.exists(path2tfa_analysis)){
        dir.create(path2tfa_analysis)
}
saveRDS(auc.df, 
        file = paste0(path2tfa_analysis, 
                     '/auc_regulons_filtered_by_size.rds'),
        compress = TRUE)


##-----------------------------------------------
## Adding information to the seurat object
all(colnames(seurat) == rownames(auc.df) )
## [1] TRUE
#seurat@meta.data <- cbind(seurat@meta.data, auc.df) 



##----------------------------------------------
## Performing differential activity analysis
tfa <- CreateAssayObject(data = auc.mtx, min.cells = 0, min.features = 0)
seurat[["regulon_"]] <- tfa
Idents(seurat) <- seurat$orig.ident
diff.tfa.list <- lapply(conditions, 
                   function(cond){
                           df <- FindMarkers(seurat,
                                       ident.1 = cond,
                                       ident.2 = '7_3h_-IFN_polyC',
                                       logfc.threshold = 0,
                                       min.pct = 0,
                                       assay = 'regulon_')
                           df$'cluster' <- cond
                           return(df)
                   })
diff.tfa.list <- diff.tfa.list[conditions != '7_3h_-IFN_polyC']
diff.tfa <- do.call(rbind, diff.tfa.list)
diff.tfa %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC) -> top



##---------------------------------------------
## Visualisation of DEA by sample
diff.tfa$'gene' <- gsub('.*\\.', '', rownames(diff.tfa))
## adding regulon size
diff.tfa$'regulon_size' <- plyr::mapvalues(diff.tfa$gene,
                                           from = reg.sizes.df$tf,
                                           to = reg.sizes.df$regulon_size)
diff.tfa$regulon_size <- as.numeric(diff.tfa$regulon_size)
## adding mean expression
gene.features <- unique(diff.tfa$gene)
seu <- NormalizeData(seurat) %>%
        ScaleData()
mean.gex.list <- lapply(conditions,
                        function(cond){
                                mtx <- FetchData(subset(seu,
                                                        orig.ident == cond),
                                                 vars = gene.features,
                                                 assay = 'RNA')
                                mean <- apply(mtx, 2, mean) 
                                df <- data.frame(gene=names(mean),
                                                 mean=mean,
                                                 condition=cond)
                                return(df)
                                
                        })
mean.gex <- do.call(rbind, mean.gex.list)
diff.tfa$'mean' <- plyr::mapvalues(rownames(diff.tfa), 
                                 from = rownames(mean.gex),
                                 to = mean.gex$mean)
diff.tfa$mean <- as.numeric(diff.tfa$mean)
vulcano.tfa <- lapply(conditions[conditions != '7_3h_-IFN_polyC'], 
                      function(cond){
        df <- filter(diff.tfa, cluster == cond)
        df <- arrange(df, desc(avg_log2FC))
        top <- head(df)$gene
        bottom <- tail(df)$gene
        df <- df %>%
                mutate(label=ifelse(gene %in% c(top, bottom,
                                                'STAT1',
                                                'STAT2'),
                                    gene, ''))
        df %>%
                ggplot(aes(x=avg_log2FC, y=-log10(p_val),
                           colour=mean,
                           label=label,
                           size=regulon_size)) +
                geom_point() +
                geom_point(data = df,
                           aes(x=avg_log2FC, 
                                    y=-log10(p_val),
                                    size=regulon_size),
                           shape=1, 
                           colour='black', 
                           alpha=0.5) +
                geom_text_repel(colour='black',
                                max.overlaps = 1000,
                                force = 10,
                                size = 4.5) +
                xlim(-0.1, 0.2) +
                scale_colour_viridis() +
                theme_classic() +
                theme(axis.text = element_text(size = 12),
                      axis.title = element_text(size = 14)) +
                labs(x='Mean Log2 FC', 
                     y='-Log10(p-value)',
                     title = substr(cond, 3, nchar(cond)))
})
pdf(paste0(path2figures, '/vulcano_plots_tfa.pdf'),
    width = 20, height = 8)
gridExtra::grid.arrange(grobs=vulcano.tfa, ncol=2)
dev.off()





##---------------------------------------------
## Plotting regulons
ntop <- 2000
regulon <- 'STAT2'


plot_rewiring <- function(regulon='STAT2', ntop=2000){
        pdf(file = paste0(path2figures, '/regulon_rewiring_', regulon, '.pdf'),
            width = 9, height = 7.5)
        par(mfrow=c(2,3))
        for ( sample in conditions){
                net.df <- networks.list[sample][[1]] %>%
                        arrange(desc(abs(coef_abs))) %>%
                        head(ntop) %>%
                        filter(source==regulon)
                edges <- select(net.df,
                                source,
                                target,
                                coef_abs,
                                coef_mean) 
                edges <- rename(edges, 
                                from=source,
                                to=target)
                nodes <- data.frame(name=unique(c(edges$from,
                                                  edges$to)))
                g <- graph_from_data_frame(edges, 
                                           directed = TRUE, 
                                           vertices = nodes)
                ## color label assignation
                nodes.label.color <- ifelse(V(g)$name == regulon, 'red', 'black')
                ## color edge assignation
                edges.color <- ifelse(E(g)$coef_mean>0, 'chartreuse3', 'orange3')
                
                ## label size
                nodes.label.size <- ifelse(V(g)$name == regulon, 1.5, 1)
                plot(g,
                     ## vertex settings
                     vertex.size=0,
                     vertex.label.color=nodes.label.color,
                     vertex.label.dist=2,
                     vertex.label.cex=nodes.label.size,
                     vertex.shape='none',
                     ## 
                     edge.arrow.size=0,
                     edge.width=(edges$coef_mean - min(edges$coef_mean, na.rm = TRUE))*10,
                     edge.color=edges.color,
                     
                     ## general settings
                     main = substr(sample, 3, nchar(sample)),
                )  
        } 
        dev.off()
}


DotPlot(seurat, features = top30.targets) +
        scale_color_viridis() +
        labs(x='', y='')

tfs <- c('IRF2', 'IRF9', 'NFKB1', 'STAT1', 'STAT2', 'JUN', 'JUNB')
for ( tf in tfs ){
        plot_rewiring(regulon = tf)
}



##---------------------------------------------
## Visualisation of dot plots
plot_dotplots <- function(regulon='STAT2', ntop=2000, genesByCondition=5){
        top.targets <- lapply(conditions, function(cond){
                top.targets <- networks.list[cond][[1]] %>%
                        arrange(desc(abs(coef_abs))) %>%
                        head(ntop) %>%
                        filter(source==regulon) %>%
                        head(genesByCondition) %>%
                        pull(target) %>%
                        unique()
        }) %>% unlist() %>% unique()
        
        
        DotPlot(seurat, features = c(paste0('sct_', top.targets)),
                assay = 'SCT') +
                scale_color_viridis() +
                labs(x='', y='', title = regulon)
}

dot.plots <- lapply(c('IRF9', 'STAT1', 'NFKB1'), 
                    function(tf) {plot_dotplots(regulon=tf, genesByCondition = 3)})
gridExtra::grid.arrange(grobs=dot.plots, ncol=1)
dot.plots[1]

##---------------------------------------------
## Visualization of TFA activities

## annotations
column_ann <- columnAnnotation(
        treatment = seurat$orig.ident,
        col = list(treatment=sample.colors[unique(seurat$orig.ident)])
)
Heatmap(
        auc.mtx[top$gene,], 
        top_annotation = column_ann, 
        row_dend_reorder = NA, 
        column_dend_reorder = TRUE,
        show_column_names = FALSE
)
dev.off()

DoHeatmap(seurat, features = top10$gene)

diff.tfa %>%
        ggplot(aes(x=avg_log2FC, y=-log10(p_val))) +
                geom_point()
