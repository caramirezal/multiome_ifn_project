## This script contain a transcription factor analysis of cells primed with IFN
## and stimulated with dsRNA

library(Seurat)
library(dplyr)
library(readr)
library(AUCell)
library(ComplexHeatmap)

## initial settings
source('/media/ag-cherrmann/cramirez/multiome_ifn_project/scripts/settings.R')

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
                geom_histogram(bins = 100) +
                geom_vline(xintercept = 30, linetype='dashed', colour='red') +
                theme_bw()



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
seurat@meta.data <- cbind(seurat@meta.data, auc.df) 



##----------------------------------------------
## Performing differential activity analysis
tfa <- CreateAssayObject(data = auc.mtx, min.cells = 0, min.features = 0)
seurat[["regulon_"]] <- tfa
Idents(seurat) <- seurat$orig.ident
diff.tfa <- FindAllMarkers(seurat, 
                        assay = 'regulon_', 
                        group.by = 'orig.ident', 
                        logfc.threshold = 0, 
                        min.pct = 0)
diff.tfa %>%
        group_by(cluster) %>%
        top_n(n = 5, wt = avg_log2FC) -> top10

seurat@assays$regulon_@scale.data <- auc.mtx



##---------------------------------------------
## Visualization of TFA activities

## annotations
column_ann <- columnAnnotation(
        treatment = seurat$orig.ident,
        col = list(treatment=sample.colors[unique(seurat$orig.ident)])
)
Heatmap(
        auc.mtx[top10$gene,], 
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
