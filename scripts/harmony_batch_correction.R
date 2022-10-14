## In this script we provide an evaluation of ATAC-Seq data after batch correction
## using harmony

## dependencies
library(ArchR)

## initial settings
set.seed(333)


## loading additional settings
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')



##-----------------------
## Loading the project
## **The archerFolder object is provided in settings.R
proj <- loadArchRProject(archerFolder)

## nthreads is provided in settings.R
addArchRThreads(threads = nthreads) 




## LSI normalization
proj <- addIterativeLSI(ArchRProj = proj, 
    useMatrix = "TileMatrix", 
    name = "Harmony",
    force = TRUE
)



##------------------------
## Correcting sample effects
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)



##-----------------------
## UMAP projection
proj <- addUMAP(ArchRProj = proj, 
        reducedDims = "Harmony",
        name = "UMAPHarmony", 
        nNeighbors = 300,
        minDist = 0.01,
        threads = nthreads,
        force = TRUE,
)


##--------------------------
## Visualization of samples
p1 <- plotEmbedding(ArchRProj = proj, 
        colorBy='cellColData',
        name = "Sample", 
        embedding = "UMAPHarmony",
        size=1.3,
        labelMeans=FALSE,
        pal=sample.colors,
)
pdf(paste0(path2project, 'figures/umap_after_harmony_batch_correction.pdf'))
p1
dev.off()






##------------------------
## Evaluation of bias due to sequencing
cells.col.df <- getCellColData(proj)
umap.df <- getEmbedding(proj)
colnames(umap.df) <- c('umap1', 'umap2')
## checking barcodes order
all(rownames(cells.col.df)==rownames(umap.df))
umap.df <- cbind(umap.df, cells.col.df)

p1 <- umap.df %>%
        ggplot(aes(x=umap1,
                   y=umap2,
                   colour=PromoterRatio),
                   size=2.5) +
                   geom_point() +
                        theme_void() +
                        scale_color_viridis() +
                        labs(title='Promoter Ratio')
p2 <- umap.df %>%
        ggplot(aes(x=umap1,
                   y=umap2,
                   colour=nFrags),
                   size=2.5) +
                   geom_point() +
                        theme_void() +
                        scale_color_viridis() +
                        labs(title='n Fragments')

pdf(paste0(figures_folder, '/umap_projection_sequencing_bias_after_harmony_batch_correction.pdf'),
        height = 6,
        width = 12
)
p1 + p2
dev.off()