## This scripts contains the annotation of the peaks in the ATAC-Seq data
## from cell lines stimulated with IFN

## dependencies
library(annotatr)
library(GenomicRanges)
library(readr)


## path2analysis
path2analysis <- paste0(path2project, 'analysis/differential_accessibility_analysis')
if ( ! dir.exists(path2analysis)){
        dir.create(path2analysis)
}

## path to annotations
path2annotations <- paste0(path2project, 'analysis/peaks_annotation')
if ( ! dir.exists(path2annotations)){
        dir.create(path2annotations)
}

peaks <- read_tsv(file = paste0(path2analysis, 
                                'differential_peaks_polyC.tsv.gz'))


dm_regions <- data.frame(seqnames=gsub(':.*', '', peaks$peak),
                         start=gsub('.*:|-.*', '', peaks$peak),
                         end=gsub('.*-', '', peaks$peak),
                         strand='*')
dm_regions <- makeGRangesFromDataFrame(dm_regions)

annots = c('hg38_basicgenes')

# Build the annotations (a single GRanges object)
annotations = build_annotations(genome = 'hg38', annotations = annots)

# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
        regions = dm_regions,
        annotations = annotations,
        ignore.strand = TRUE,
        quiet = FALSE)
# A GRanges object is returned
print(dm_annotated)


# Coerce to a data.frame
df_dm_annotated = data.frame(dm_annotated)

# See the GRanges column of dm_annotaed expanded
print(head(df_dm_annotated))
dim(df_dm_annotated)

annotated_unambigously <- df_dm_annotated[
        !duplicated(dplyr::select(df_dm_annotated,
                                  seqnames,
                                  start,
                                  end)),
]


head(annotated_unambigously)
dim(annotated_unambigously)


saveRDS(annotated_unambigously,
        file = paste0(path2annotations, '/annotated_peaks_unambiguously.rds'), 
        compress = TRUE)
