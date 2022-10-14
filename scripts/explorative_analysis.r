## In this script we provide an explorative data analysis of single cell multiomic data
## from human cells treated or not with IFN and polyC

## dependencies
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(viridis)

##--------------------------
## Initial settings
path2project <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/'
setwd(path2project <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/')
set.seed(333)

## loading additional general settings 
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')

## threads
nthreads <- 32
addArchRThreads(threads = nthreads) 

## re-run or load previous results
run_archer_project_pipeline <- FALSE

sample.colors <- c(
        '7_3h_-IFN_polyC' = 'chartreuse',
        '6_3h_-IFN_dsRNA' = 'chartreuse3',
        '5_3h_pIFN_polyC' = 'lightsalmon',
        '4_3h_pIFN_dsRNA' = 'darkorange',
        '1_16h_pIFN_dsRNA'= 'firebrick3',
        '2_16h_pIFN_polyC'= 'firebrick4',
        '3_16h_-IFN_dsRNA'= 'chartreuse4' 
)

##--------------------------
## Setting inputs
fragment.files <- list.files(
        '/media/sds-hd/sd21e005/binder_multiome/counts/', 
        pattern = 'atac_fragments.tsv.gz$', 
        full.names = TRUE, 
        recursive = TRUE
)
names(fragment.files) <- gsub('.*//', '', fragment.files)
names(fragment.files) <- gsub('/.*', '', names(fragment.files)) 
fragment.files


##################################################################################
##                                                                              ##
##                      Archer project pipeline                                 ##
##`                                                                             ##
##################################################################################

if ( run_archer_project_pipeline == TRUE ) {

## Setting genome reference
addArchRGenome("hg38")

##-------------------------
## Creation of arrow files

ArrowFiles <- createArrowFiles(
        inputFiles = fragment.files,
        sampleNames = names(fragment.files),
        filterTSS = 4, #Dont set this too high because you can always increase later
        filterFrags = 1000, 
        addTileMat = TRUE,
        addGeneScoreMat = TRUE,
)

##--------------------------
## Inference of doublets
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
)


##------------------------
## Creating archer project
archerFolder <- '/media/sds-hd/sd21e005/binder_multiome/archer_folder'
dir.create(archerFolder)
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = archerFolder,
  copyArrows = TRUE, #This is recommened so that you maintain an unaltered copy for later usage.
)
getAvailableMatrices(proj)



##-----------------------
## filtering doublets
proj <- filterDoublets(ArchRProj = proj)



##-----------------------
## Dimensional reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")



##-----------------------
## Clustering
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")



##-----------------------
## UMAP projection
proj <- addUMAP(ArchRProj = proj, 
        reducedDims = "IterativeLSI", 
        nNeighbors = 50,
        minDist = 0.8,
        threads = nthreads,
        force = TRUE,
)

## Saving the project
saveArchRProject(
        ArchRProj = proj,
        outputDirectory = archerFolder,
        overwrite = TRUE,  
)



}  ## End of the Archer project pipeline


######################################################################

loadArchRProject(path=archerFolder, force=FALSE, showLogo = TRUE)

##---------------------
## Plotting UMAPs

## Visualization of samples
p1 <- plotEmbedding(ArchRProj = proj, 
        colorBy='cellColData',
        name = "Sample", 
        embedding = "UMAP",
        size=1.3,
        labelMeans=FALSE,
        pal=sample.colors,
)

figures_folder <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/'
dir.create(figures_folder)
pdf(paste0(figures_folder, '/umap_projection_all_samples.pdf'),
        height = 6,
        width = 6,
)
p1 + theme_classic()
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

pdf(paste0(figures_folder, '/umap_projection_sequencing_bias.pdf'),
        height = 6,
        width = 12
)
p1 + p2
dev.off()




##-------------------------------------
## Visualization of IFN markers
proj <- addImputeWeights(proj)

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = 'ISG15', 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj),
    pal = viridis()
)

pdf(paste0(figures_folder, '/umap_isg15_expression.pdf'), 
        height = 5,
        width = 5
)
p
dev.off()


