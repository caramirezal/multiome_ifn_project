## Peak calling by sample


## dependencies
library(ArchR)
library(ggplot2)
library(viridis)
library(Matrix)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(GenomicRanges)
library(pheatmap)

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
run_peak_calling <- TRUE


proj <- loadArchRProject(path=archerFolder, force=FALSE, showLogo = TRUE)


## Peak calling by sample
samples <- unique(proj$Sample)
names(samples) <- samples
if ( run_peak_calling == TRUE ) {
        proj.samp <- lapply(samples, 
                            function(samp){
                                    
                                    proj.samp <- proj[proj$Sample == samp, ]
                                    
                                    ##-----------------------
                                    ## Dimensional reduction
                                    proj.samp <- addIterativeLSI(
                                            ArchRProj = proj.samp, 
                                            useMatrix = "TileMatrix", 
                                            name = "IterativeLSI", 
                                            force = TRUE
                                    )
                                    
                                    
                                    
                                    ##-----------------------
                                    ## Clustering
                                    proj.samp <- addClusters(
                                            input = proj.samp, 
                                            reducedDims = "IterativeLSI",
                                            force=TRUE
                                    )
                                    
                                    
                                    
                                    ##-----------------------
                                    ## Creation of pseudo-bulk replicates as controls for
                                    ## the peak calling
                                    proj.samp <- addGroupCoverages(
                                            ArchRProj = proj.samp, 
                                            groupBy = "Clusters",
                                            force = TRUE
                                    )
                                    
                                    
                                    
                                    ## Calling peaks in ATAC Seq data
                                    pathToMacs2 <- findMacs2()
                                    proj.samp <- addReproduciblePeakSet(
                                            ArchRProj = proj.samp, 
                                            groupBy = "Clusters", 
                                            pathToMacs2 = pathToMacs2
                                    )
                                    proj.samp <- addPeakMatrix(proj.samp)
                                    
                                    
                                    ## Saving results
                                    saveRDS(proj.samp, '/media/sds-hd/sd21e005/binder_multiome/archer_folder/peaksBysample.rds')
                                    
                                    
                            })
        
}



##--------------------
## Reading pre-processed peaks by sample
proj.samp <- readRDS('/media/sds-hd/sd21e005/binder_multiome/archer_folder/peaksBysample.rds')




##------------------------
## Extracting the peaks
ranges.df <- lapply(proj.samp, function(proj){
        getMatrixFromProject(proj, useMatrix = 'PeakMatrix') %>% 
                rowRanges()
})




##------------------------
## Definition of a matrix to store the number of peak overlaps
overlap.count.matrix <- matrix(data = 0, 
                               nrow = length(samples),
                               ncol = length(samples))
colnames(overlap.count.matrix) <- samples
rownames(overlap.count.matrix) <- samples



## -----------------------
## Counting overlaps
for (i in rownames(overlap.count.matrix)){
        for ( j in colnames(overlap.count.matrix)) {
                overlap.count.matrix[i, j] <- countOverlaps(query = ranges.df[i][[1]], 
                                                             subject = ranges.df[j][[1]]) %>%
                                                        sum()
        }
}
        


##--------------------------
## Dividing number of overlaps by the total number of peaks
overlap.count.matrix.norm <- overlap.count.matrix / diag(overlap.count.matrix)
#overlap.count.matrix.norm <- t(scale(overlap.count.matrix.norm))
#diag(overlap.count.matrix.norm) <- 0



## Printing heatmap with the number of overlaps
pdf('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/heatmap_overlaps.pdf',
    height = 10, width = 10)
pheatmap(overlap.count.matrix.norm, cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()



##------------------------------
## Counting total number of peaks
npeaks <- sapply(ranges.df, function(grange) nrow(as.data.frame(grange)))
npeaks.df <- data.frame(sample = names(npeaks), nb_peaks=npeaks)



## barplot of the number of peaks
pdf('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/barplot_no_peaks.pdf')
npeaks.df %>%
        arrange(desc(nb_peaks)) %>%
        mutate(sample=factor(sample, levels = sample)) %>%
        ggplot(aes(x=sample, y= nb_peaks)) +
                geom_bar(stat = 'identity', fill='steelblue') +
                theme_bw() +
                theme(panel.grid = element_blank(),
                      axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(x='', y='Number of peaks') 
dev.off()

