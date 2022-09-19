## This script contains the bioinformatic pipeline used to process multiomic data
## from human cell lines Treated/Untreated with IFN or polyIC

path2project <- paste0('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/')

## First the preprocessing_fastq.R. Fastq files need to be renamed in order to fulfill
## cellranger requirements. With this scripts hardlinks are created and additionally
## the libraries files necessary to define the cell ranger counting are defined
## Inputs:
## - ATACSeq metadata: xlsx file in data/30773-resultData.xls
## - RNASeq metadata: tsv file in data/data/220812_A00382_0469_AHHWGMDRX2/220812_A00382_0469_AHHWGMDRX2_meta.tsv
## Outputs:
## - libraries: stored in data/libraries
## - hardlinks: stored in /media/sds-hd/sd21e005/binder_multiome/data/hardlinks
source(paste0(path2project, '/scripts/preprocessing_fastq.R'))


## A check on the fastq integrity is given in the checksum_check.R file. This
## step is optional. A table containing the checksums of the cell ranger output
## and the downloaded files are provided in the table
## outputs: 
## - checksum table: data/checksums/checksums.xlsx
source(paste0('/scripts/checksum_check.R'))


## run_counting.sh counts the reads performing the cellranger-arc count algorithm.
## This file can be executed interactively or using the qsub_job.sh in the
## qsub cluster.
## Inputs:
## - library files from above
## - fastq files (hard links) defined above
## Outputs:
## - counts: stored in '/media/sds-hd/sd21e005/binder_multiome/counts'
system(paste0(path2project, 'scripts/qsub_job.sh'))


