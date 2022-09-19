## This script contains the pre-processing and counting fastq files of RNA and ATAC multiomic data
## of cell cultures treated or not with IFN and/or polyI:C

## dependencies
library(tidyverse)

##-----------------------
## Preprocessing RNA Seq data

## Reading RNA Seq metadata
path2rna_metadata <- '/media/sds-hd/sd21e005/binder_multiome/data/data/220812_A00382_0469_AHHWGMDRX2/220812_A00382_0469_AHHWGMDRX2_meta.tsv'
rna_meta <- read_tsv(path2rna_metadata)

rna_meta %>% as.data.frame() 

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

## Defining new folder for the hardlinks
path2rna_new <- '/media/sds-hd/sd21e005/binder_multiome/data/hardlinks/rna/'
## Creating new folders
folders.old <- list.dirs('/media/sds-hd/sd21e005/binder_multiome/data/data/220812_A00382_0469_AHHWGMDRX2',
          recursive = FALSE)
folders.old <- gsub('.*/', '', folders.old)
folders.old <- folders.old[!folders.old=="renamed_rna"]
for ( folder in paste0(path2rna_new, folders.old) ){
        if ( ! dir.exists(folder) ) {
                dir.create(folder, recursive = TRUE)
        }
}
list.dirs(path2rna_new)

## renaming files
rna_meta <- mutate(rna_meta, fastq_path_renamed=paste0(path2rna_new, 
                                                       gsub('_.*', '', FASTQ_FILE))) 

## old paths
path2rna <- '/media/sds-hd/sd21e005/binder_multiome/data/data/220812_A00382_0469_AHHWGMDRX2/'
rna_meta <- mutate(rna_meta, fastq_path=paste0(path2rna, 
                                                       gsub('_.*', '', FASTQ_FILE), '/',
                                                       'fastq/',
                                               FASTQ_FILE))
rna_meta

##------------------
## Saving renamed fastq files
write_tsv(rna_meta, file = '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/rna_renamed_fastq.tsv')

##------------------
## Creating hard links
for (i in 1:nrow(rna_meta)) {
        command <- paste0('ln ', rna_meta$fastq_path[i], 
                          ' ', 
                          paste0(rna_meta$fastq_path_renamed[i], '/',
                                 gsub('\\+', 'p', rna_meta$fastq_renamed[i])
                                 )
                          ) 
        cat('Executing', command, '\n')
        system(command)
}
list.files(path2rna_new, recursive = TRUE)




##----------------------
## Preprocessing ATAC Seq data

## Reading ATAC Seq data
path2atac <- '/media/sds-hd/sd21e005/binder_multiome/data/data/220805_A00382_0466_AHFCYTDMXY/'
atac.fastq.files <- list.files(path2atac, pattern = 'fastq.gz$', recursive = TRUE)
atac.fastq.files <- atac.fastq.files[!grepl('Undetermined|mySample', atac.fastq.files)]
atac.fastq.files <- atac.fastq.files[!grepl('_S1_', atac.fastq.files)]

##---------------------
## Reading ATAC seq metadata
atac_metadata_path <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/30773-resultData.xls'
atac_info <- readxl::read_xls(atac_metadata_path)
atac_info <- rename(atac_info, id_lane = `Unique ID / Lane`)

## Adding sample information
atac_meta <- data.frame(fastq_file=atac.fastq.files)
atac_meta <- mutate(atac_meta, id_lane=gsub('/.*', '', fastq_file))
atac_meta <- merge(atac_meta, atac_info)

head(atac_meta)
atac_meta <- select(atac_meta, fastq_file, `Sample Name`)
atac_meta <- mutate(atac_meta, fastq_file_path = paste0(path2atac, fastq_file),
                               sample_name = gsub('_[A-D]$', '', `Sample Name`))
## checking that file exists
all(sapply(atac_meta$fastq_file_path, file.exists))

## lane
atac_meta <- mutate(atac_meta, lane_number=gsub('.*-LR-|_[IR][12].*', '', fastq_file))
atac_meta <- mutate(atac_meta, lane_number=as.integer(as.factor(lane_number)))

## read type
atac_meta <- mutate(atac_meta, read_type=gsub('.fastq.gz', '', fastq_file))
atac_meta <- mutate(atac_meta, read_type=gsub('.*_', '', read_type))

## sample id + lane
atac_meta <- mutate(atac_meta, sample_lane=gsub('/fastq.*', '', fastq_file))

        
## Renaming
atac_meta <- mutate(atac_meta,
                    fastq_renamed=paste0(
                            sample_name,
                            '_S1_L00',
                            lane_number,
                            '_', read_type,
                            '_001.fastq.gz'
                    ))
atac_meta <- mutate(atac_meta, fastq_renamed=gsub('\\+', 'p', fastq_renamed))

## Defining new paths
path2atac_new <- '/media/sds-hd/sd21e005/binder_multiome/data/hardlinks/atac/'
atac_meta <- mutate(atac_meta, 
                    fastq_file_renamed_path = paste0(path2atac_new, 
                                                     sample_lane))
head(atac_meta)

write_tsv(atac_meta, file = '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/atac_renamed_fastq.tsv')

## Creating folders for the atac hard links
atac_folders <- atac_meta$sample_lane %>% unique()
atac_folders <- paste0(path2atac_new, atac_folders)
for (folder in atac_folders) {
        if ( ! dir.exists(folder)){
                dir.create(folder)
        }
}

##------------------
## Creating hard links
for (i in 1:nrow(atac_meta)) { 
        command <- paste0('ln ', atac_meta$fastq_file_path[i], ' ', 
                          paste0(atac_meta$fastq_file_renamed_path[i],'/', atac_meta$fastq_renamed[i] ) )
        cat('Executing', command, '\n')
        #system(command)
}

list.files(path2atac_new, pattern = 'fastq.gz$', recursive = TRUE)



##--------------------------
## Merging information


rna_libraries <- select(rna_meta, fastq_path_renamed, SAMPLE_NAME) 
rna_libraries <- mutate(rna_libraries, library_type='Gene Expression')
rna_libraries <- rename(rna_libraries, fastqs=fastq_path_renamed, sample=SAMPLE_NAME)

atac_libraries <- select(atac_meta, fastq_file_renamed_path, sample_name)
atac_libraries <- mutate(atac_libraries, library_type='Chromatin Accessibility')
atac_libraries <- rename(atac_libraries, sample=sample_name, fastqs=fastq_file_renamed_path)

multi_libraries <- rbind(rna_libraries, atac_libraries)
multi_libraries <- mutate(multi_libraries, sample=gsub('\\+', 'p', sample))
samples <- unique(multi_libraries$sample)

path2librariesCSV <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/libraries'
if ( ! dir.exists(path2librariesCSV)){
        dir.create(path2librariesCSV)
}


for (samp in samples) {
        multi_library.df <- filter(multi_libraries, sample==samp)
        multi_library.df <- unique(multi_library.df)
        cat('Processing', paste0(path2librariesCSV, '/', samp, '_libraries.csv'), '\n')
        write_csv(multi_library.df, 
                  file = paste0(path2librariesCSV, '/', samp, '_libraries.csv'))
}





