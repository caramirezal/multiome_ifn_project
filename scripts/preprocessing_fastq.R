## This script contains the pre-processing and counting fastq files of RNA and ATAC multiomic data
## of cell cultures treated or not with IFN and/or polyI:C

library(tidyverse)

##-----------------------
## Preprocessing RNA Seq data

## Reading RNA Seq metadata
path2rna_metadata <- '/media/sds-hd/sd21e005/binder_multiome/data/data/220812_A00382_0469_AHHWGMDRX2/220812_A00382_0469_AHHWGMDRX2_meta.tsv'
rna_meta <- read_tsv(path2rna_metadata)

## Adding hard link name
rna_meta <- filter(rna_meta, !grepl('Undetermined', FASTQ_FILE))
rna_meta <- mutate(rna_meta, 
       fastq_renamed = paste0(SAMPLE_NAME,
                              '_S1_L00',
                              LANE_NO,
                              '_',
                              'R', READ,
                              '_001.fastq.gz'
               
       ))
path2rna <- '/media/sds-hd/sd21e005/binder_multiome/data/data/220812_A00382_0469_AHHWGMDRX2/'
rna_meta <- mutate(rna_meta, fastq_path_renamed=paste0(path2rna, 
                                                       gsub('_R.*', '', FASTQ_FILE), '/fastq/',
                                                       fastq_renamed)) 
rna_meta <- mutate(rna_meta, fastq_path=paste0(path2rna, 
                                                       gsub('_R.*', '', FASTQ_FILE), '/fastq/',
                                                       FASTQ_FILE)) 

## Creating hard links
for (i in 1:nrow(rna_meta)) { 
        command <- paste0('ln ', rna_meta$fastq_path[i], ' ', rna_meta$fastq_path_renamed[i]) 
        cat('Executing', command, '\n')
        system(command)
}


