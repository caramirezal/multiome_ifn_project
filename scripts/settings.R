## In this file we provide settings for the overal project, for example,
## categories colors, ...

## project github path
path2project <- '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/'

## Path to archer project
archerFolder <- '/media/sds-hd/sd21e005/binder_multiome/archer_folder'


## sample colors
sample.colors <- c(
        '7_3h_-IFN_polyC' = 'chartreuse',
        '6_3h_-IFN_dsRNA' = 'chartreuse3',
        '5_3h_pIFN_polyC' = 'lightsalmon',
        '4_3h_pIFN_dsRNA' = 'darkorange',
        '1_16h_pIFN_dsRNA'= 'firebrick3',
        '2_16h_pIFN_polyC'= 'firebrick4',
        '3_16h_-IFN_dsRNA'= 'chartreuse4' 
)

## Paralellization
nthreads <- 9