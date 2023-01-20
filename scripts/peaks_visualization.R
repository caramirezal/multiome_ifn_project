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
pdf(paste0(figures_folder, '/peaksTracks_upregulated_after_priming.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '7_3h_-IFN_polyC'), ], 
        groupBy = "Sample", 
        geneSymbol = c("REEP6", "SAMD4B", "CIC"),
        features =  getPeakSet(proj),
        upstream = 30000,
        downstream = 30000,
        pal=sample.colors
)
grid::grid.draw(p$REEP6)
grid::grid.newpage()
grid::grid.draw(p$SAMD4B)
grid::grid.newpage()
grid::grid.draw(p$CIC)
dev.off()




print('Plotting')
pdf(paste0(figures_folder, '/peaksTracks_downregulated_after_priming.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '7_3h_-IFN_polyC'), ], 
        groupBy = "Sample", 
        geneSymbol = c("PLEKHH2", "GRK6", "ZNRD2"),
        features =  getPeakSet(proj),
        upstream = 30000,
        downstream = 30000,
        pal=sample.colors
)
grid::grid.draw(p$PLEKHH2)
grid::grid.newpage()
grid::grid.draw(p$GRK6)
grid::grid.newpage()
grid::grid.draw(p$ZNRD2)
dev.off()




print('Plotting')
pdf(paste0(figures_folder, '/peaksTracks_upregulated_after_stimulation.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '4_3h_pIFN_dsRNA'), ], 
        groupBy = "Sample", 
        geneSymbol = c("ZNF837", 'ZNF84', 'MYH10'),
        features =  getPeakSet(proj),
        upstream = 30000,
        downstream = 30000,
        pal=sample.colors
)
grid::grid.draw(p$ZNF837)
grid::grid.newpage()
grid::grid.draw(p$ERCC1)
grid::grid.newpage()
grid::grid.draw(p$ZNF84)
dev.off()


print('Plotting')
pdf(paste0(figures_folder, '/peaksTracks_downregulated_after_stimulation.pdf'), 
    width = 6, height = 6)
p <- plotBrowserTrack(
        ArchRProj = proj[proj$Sample %in% c('5_3h_pIFN_polyC', '4_3h_pIFN_dsRNA'), ], 
        groupBy = "Sample", 
        geneSymbol = c("NR1H2", 'CIC', 'PNKD'),
        features =  getPeakSet(proj),
        upstream = 30000,
        downstream = 30000,
        pal=sample.colors
)
grid::grid.draw(p$NR1H2)
grid::grid.newpage()
grid::grid.draw(p$CIC)
grid::grid.newpage()
grid::grid.draw(p$PNKD)
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