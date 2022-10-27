## This script contains the bioinformatic pipeline used to process multiomic data
## from human cell lines Treated/Untreated with IFN or polyIC

## path to github project path in the curry cluster
path2project <- paste0('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/')


## Fastq processing
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


## explorative_analysis.r contains processing of ATAC-Seq data using ArchR pipeline
## Inputs:
## - fragment files in /media/sds-hd/sd21e005/binder_multiome/counts 
## Outputs:
## - ArrowFiles stored in /media/sds-hd/sd21e005/binder_multiome/archer_folder
## - ArrowProject output in the same folder above
## - Figures: umap_projection_all_samples.pdf, umap_projection_sequencing_bias.pdf, umap_isg15_expression.pdf
system(paste0(path2project, 'scripts/explorative_analysis.R'))

## The script scripts/explorative_analysis_filtering_conditions.r provides
## an alternative analysis subsetiing to samples with good appeareance in the QC

## harmony_batch_correction.R contains an alternative processing with the ArchR pipeline
## in which harmony is used to correct for sample effects
## Inputs:
## - fragment files in /media/sds-hd/sd21e005/binder_multiome/counts 
## Outputs:
## - Figures: umap_after_harmony_batch_correction.pdf, umap_projection_sequencing_bias_after_correction.pdf
system(paste0(path2project, 'scripts/harmony_batch_correction.R'))



## rna_seq_processed_individually.R performs a analysis of only the
## gene expression data separately
## Inouts: 
## - Counts matrix: stored in the counts_GEX folder
## Output:
## - Seurat object: containing gene expression data along with annotations
## - plots: vln_plots_QC, umap_conditions_n_qc_metrics_rnaseq.pdf
system(paste0(path2project, 'scripts/rna_seq_processed_individually.R '))


## ifn_signatures evaluates the expression of signatures associated to the
## IFN response
## Inputs:
## - Seurat object: processed using the rna_seq_processed_individually.R above
## Outputs:
## - plots:
## - DEGs:


