## This script check the integrity of fastq files by comparing the checksum files

## Reading calculated check sums
md5sums <- read.table('/media/sds-hd/sd21e005/binder_multiome/md5sum_checks.txt')
colnames(md5sums) <- c('md5sum', 'fastq')
md5sums <- mutate(md5sums, fastq_path=fastq)
md5sums <- mutate(md5sums, fastq=gsub('.*/fastq/', '', fastq))


md5sums10x <- read.table('/media/sds-hd/sd21e005/binder_multiome/md5sum_checks_10x.txt')
colnames(md5sums10x) <- c('md5sum_10x', 'fastq')
md5sums10x <- mutate(md5sums10x, fastq=gsub('.md5sum$', '', fastq))

head(md5sums); head(md5sums10x)
checksums <- merge(md5sums, md5sums10x)
checksums

nonPassingChecksums <- filter(checksums, md5sum != md5sum_10x)

dir.create('/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/checksums')
writexl::write_xlsx(checksums, path = '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/checksums/checksums.xlsx')
writexl::write_xlsx(nonPassingChecksums, path = '/media/sds-hd/sd21e005/binder_multiome/multiome_ifn_project/data/checksums/checksums_non-matching.xlsx')