## This script contain code to visualise ATAC tracks using ArchR

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



## Defining folder to store ouputs
figures_folder <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/figures/peaks_visualization'
if ( ! dir.exists(figures_folder)){
        dir.create(figures_folder)
}



proj <- loadArchRProject(path=archerFolder, force=FALSE, showLogo = TRUE)




#markerList <- readRDS(file = paste0(path2project, 'analysis/differential_peaks.rds'))

print('Plotting')
pdf(paste0(figures_folder, '/peaksTracks_REEP6.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '7_3h_-IFN_polyC'), ], 
        groupBy = "Sample", 
        geneSymbol = c("REEP6"),
        features =  getPeakSet(proj),
        upstream = 30000,
        downstream = 30000,
        pal=sample.colors
)
grid::grid.draw(p$REEP6)
dev.off()




print('Plotting')
pdf(paste0(figures_folder, '/peaksTracks_ZNF837.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '4_3h_pIFN_dsRNA'), ], 
        groupBy = "Sample", 
        geneSymbol = c("ZNF837"),
        features =  getPeakSet(proj),
        upstream = 30000,
        downstream = 30000,
        pal=sample.colors
)
grid::grid.draw(p$ZNF837)
dev.off()


print('Plotting')
pdf(paste0(figures_folder, '/peaksTracks_BCL9L.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '7_3h_-IFN_polyC', 
                                            '4_3h_pIFN_dsRNA'), ], 
        groupBy = "Sample", 
        geneSymbol = c("BCL9L"),
        features =  getPeakSet(proj),
        upstream = 5000,
        downstream = 5000,
        pal=sample.colors
)
grid::grid.draw(p$BCL9L)
dev.off()