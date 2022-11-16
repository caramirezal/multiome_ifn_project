## In this script we provide an explorative data analysis of single cell multiomic data
## from human cells treated or not with IFN and polyC
## In this analysis we subset the samples to that with good appeareance in the QC

## dependencies
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(grid)



##--------------------------
## Initial settings
path2project <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/'
setwd(path2project <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/')
set.seed(333)

## loading additional general settings 
source('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/scripts/settings.R')



## Folder to store archR files
archerFolder <- '/media/sds-hd/sd21e005/binder_multiome/archer_folder'
if ( ! dir.exists(archerFolder)){
        dir.create(archerFolder)
}

## Defining folder to store ouputs
figures_folder <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/explorative_analysis_all_conditions'
if ( ! dir.exists(figures_folder)){
        dir.create(figures_folder)
}



## threads
nthreads <- 32
addArchRThreads(threads = nthreads) 

## re-run or load previous results
run_archer_project_pipeline <- TRUE

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
## Subsetting to samples with good quality
## appeareance
subset <- names(sample.colors)

sample.colors <- sample.colors[names(sample.colors) %in% subset]

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


fragment.files <- fragment.files[names(fragment.files) %in% subset]

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
if ( ! dir.exists(archerFolder)){
     dir.create(archerFolder)
}
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
## Removing cells with high number of fragments > 5000
## We observed a bias based in the UMAP representation
proj <- proj[proj$nFrags < 5000, ]



##-----------------------
## Dimensional reduction
proj <- addIterativeLSI(
        ArchRProj = proj, 
        useMatrix = "TileMatrix", 
        name = "IterativeLSI", 
        force = TRUE
)



##-----------------------
## Clustering
proj <- addClusters(
        input = proj, 
        reducedDims = "IterativeLSI",
        force=TRUE
)



##-----------------------
## Creation of pseudo-bulk relpicates as controls for
## the peak calling
proj <- addGroupCoverages(
        ArchRProj = proj, 
        groupBy = "Clusters",
        force = TRUE
)



## Calling peaks in ATAC Seq data
pathToMacs2 <- findMacs2()
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)




##-----------------------
## Dimensional reduction
proj <- addIterativeLSI(
        ArchRProj = proj, 
        useMatrix = "PeakMatrix", 
        name = "IterativeLSI", 
        force = TRUE
)




##-----------------------
## UMAP projection
proj <- addUMAP(ArchRProj = proj, 
        reducedDims = "IterativeLSI", 
        nNeighbors = 200,
        minDist = 0.0001,
        threads = nthreads,
        force = TRUE,
)


##-----------------------
## Saving the project
saveArchRProject(
        ArchRProj = proj,
        outputDirectory = archerFolder,
        overwrite = TRUE,  
)



}  ## End of the Archer project pipeline


######################################################################

proj <- loadArchRProject(path=archerFolder, force=FALSE, showLogo = TRUE)





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



pdf(paste0(figures_folder, '/umap_projection_using_peaks_nfrags<5000.pdf'),
        height = 6,
        width = 6,
)
p1 + theme_classic() + labs(colour='')
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

p3 <- plotEmbedding(ArchRProj = proj, 
                    colorBy='cellColData',
                    name = "Clusters", 
                    embedding = "UMAP",
                    size=1.3,
                    labelMeans=FALSE
)


cat('Plotting UMAP to assess sequencing bias')
pdf(paste0(figures_folder, '/umap_projection_sequencing_bias_using_peaks_nfrags<5000.pdf'),
        height = 6,
        width = 6
)
p1
p2
p3
dev.off()





##########################################################################################
## Differential peak analysis



##--------------------------
## Calculating differential peaks
cat('Calculating differential peaks')
markersPeaks <- getMarkerFeatures(
        ArchRProj = proj, 
        useMatrix = "PeakMatrix", 
        groupBy = "Sample",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 1")
saveRDS(markerList, file = paste0(path2project, 'analysis/differential_peaks.rds'))


markerList <- readRDS(file = paste0(path2project, 'analysis/differential_peaks.rds'))


pdf(paste0(figures_folder, '/fdrVsLog2FC_peaks.pdf'), width = 6, height = 6)
samples <- names(sample.colors)
names(samples) <- samples
lapply(samples, 
       function(samp){
               markerList[samp][[1]] %>%
                       as.data.frame() %>%
                       ggplot(aes(x=Log2FC, y=FDR)) +
                       geom_point() +
                       theme_bw() +
                       theme(panel.background = element_blank()) +
                                labs(title = samp)
       })
dev.off()


pdf(paste0(figures_folder, '/peaksTracks.pdf'), width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj, 
        groupBy = "Sample", 
        geneSymbol = c("ISG15"),
        features =  makeGRangesFromDataFrame(markerList$'5_3h_pIFN_polyC'),
        upstream = 50000,
        downstream = 50000,
        pal=sample.colors
)
grid::grid.draw(p$ISG15)
dev.off()

#makeGRangesFromDataFrame(markerList$'5_3h_pIFN_polyC')
