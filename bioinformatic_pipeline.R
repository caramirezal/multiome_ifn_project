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
source(paste0(path2project, '/scripts/checksum_check.R'))


## run_counting.sh counts the reads performing the cellranger-arc count algorithm.
## This file can be executed interactively or using the qsub_job.sh in the
## qsub cluster.
## Inputs:
## - library files from above
## - fastq files (hard links) defined above
## Outputs:
## - counts: stored in '/media/sds-hd/sd21e005/binder_multiome/counts'
system(paste0(path2project, '/scripts/qsub_job.sh'))


## explorative_analysis.r contains processing of ATAC-Seq data using ArchR pipeline
## in order to run this script is necessary first to build the environment in
## envs/multiome_ifn_project.yml
## Similar analysis are done removing some conditions or using all conditions and are shown in
## explorative_analysis_all_conditions.r 
## Inputs:
## - fragment files in /media/sds-hd/sd21e005/binder_multiome/counts 
## Outputs:
## - ArrowFiles stored in /media/sds-hd/sd21e005/binder_multiome/archer_folder
## - ArrowProject output in the same folder above
## - Figures: umap_projection_all_samples.pdf, umap_projection_sequencing_bias.pdf, umap_isg15_expression.pdf
system(paste0(path2project, '/scripts/explorative_analysis.R'))

## The script scripts/explorative_analysis_filtering_conditions.r provides
## an alternative analysis subsetiing to samples with good appeareance in the QC

## harmony_batch_correction.R contains an alternative processing with the ArchR pipeline
## in which harmony is used to correct for sample effects
## Inputs:
## - fragment files in /media/sds-hd/sd21e005/binder_multiome/counts 
## Outputs:
## - Figures: umap_after_harmony_batch_correction.pdf, umap_projection_sequencing_bias_after_correction.pdf
system(paste0(path2project, '/scripts/harmony_batch_correction.R'))



## rna_seq_processed_individually.R performs a analysis of only the
## gene expression data separately
## Inputs: 
## - Counts matrix: stored in the counts_GEX folder
## Output:
## - Seurat object: containing gene expression data along with annotations
## - plots: vln_plots_QC, umap_conditions_n_qc_metrics_rnaseq.pdf
system(paste0(path2project, '/scripts/rna_seq_processed_individually.R'))


## ifn_signatures evaluates the expression of signatures associated to the
## IFN response
## Inputs: 
## - Seurat object processed using the rna_seq_processed_individually.R above
## - Signatures stored in data/signatures/ifn_signatures.xlsx and the msigdb signature
##      HECKER_IFNB1_TARGETS
## Outputs:
## - AUC scores: analysis/ifn_signatures/dsRNAVspolyC_across_conditions.tsv
## - UMAP plot projection showing the IFN signature scores ('figures/inf_signatures/ifn_signatures.pdf')
## - DEGs plot figures/inf_signatures/scat_plot_comparing_IFN_signatures.pdf showing the
##    signatures comparing IFN+ dsRNA Vs IFN+ polyC and on the other hand IFN- dsRNA Vs IFN- polyC. The
##    table containing these two comparissons are stored in analysis/inf_signatures/dsRNAVspolyC_across_conditions.tsv
##   Also the signatures comparing IFN+ polyC Vs IFN- polyC and on the other hand comparing IFN+ dsRNA 
##   and IFN- dsRNA are plotted in the scat_plot_comparing_IFN_signatures_preVsNo_treatment.pdf
##   and the table are stored in the file preVsnonIFN_across_conditions.tsv
## - dot plots of the Gene Expression of genes associated to the IFN signature ('figures/ifn_signatures/ifn_signatures_by_sample.pdf',
##   figures/inf_signatures/ifn_genes_by_sample.pdf) 
system(paste0(path2project, '/scripts/ifn_signatures.R'))


## ifn_signature_corregulated_genes.R file contains an analysis on the
## corregulated genes to the IFN signatures using only the GEX data
## Inputs:
## - Seurat object processed using the rna_seq_processed_individually.R above
## - Interferon signatures list as those used in the ifn_signatures.R analysis
## Outputs:
## - ifn_signature_corregulated_genes/ifn_signatures_additional_signatures.pdf containing a umap
##  plot showing the ifn scores values
## - The output from a boosting model using as independent variables the ifn signatures scores
## and the gene expression as explanatory variables stored in the 
## ifn_signature_corregulated_genes/ifn_signature_corregulated_genes.pdf
system(paste0(path2project, '/scripts/ifn_signature_corregulated_genes.R'))


## scripts/differential_accesibility_peaks.R performs a wilcoxon test of the peak counts
## for different conditions
## Input:
## - Peak counts from the peak calling performed by ArchR and stored in the folder data/peaks
## Output:
## - Seurat object containing the peaks counts stored in the 
## analysis/differential_accessibility_analysis/atac_peaks_seu.rds
## - A DEA comparing 4_3h_pIFN_dsRNA Vs 6_3h_-IFN_dsRNA stored in 
## analysis/differential_accessibility_analysis/dea/differential_peaks_dsRNA.tsv.gz
## - A DEA comparing 5_3h_pIFN_polyC Vs 7_3h_-IFN_polyCstored in 
## analysis/differential_accessibility_analysis/dea/differential_peaks_polyC.tsv.gz
## - umap_diff_peaks.pdf A umap plot constructed by using only the Differential peaks
## calculated above
## - differential_peak_analysis_IFN(+Vs-)_dsRNa_Vs_polyIC.pdf a scatter plot with the
## AvgLog2FC comparing IFN (+Vs-) tested separately in dsRNA and polyC conditions.
system(paste0(path2project, '/scripts/differential_accesibility_analysis.R'))

