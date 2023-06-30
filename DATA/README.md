# Dataframes and input

## A549

### mapping_logs

+ Output mapping logs from alignment of sequenced reads generated from A549 ATAC-seq experiments

### mapping_stats

+ Output mapping statistics on A549 ATAC-seq experiments

## BEDS

+  Genomic counts ber 10 kb of whole fragments per kilobase per million for A549 and ENCODE samples

## ENCODE

+ logs - Output folder with logs per sample (take from ENCODE project) from alignment pipeline (described in Roth et al).

+ stats - A directory holding the calculated mapping statistics for each ENCODE sample. 

+ data.names.txt - Names of used ENCODE samples

## IGV

+ TSV files for chromosomes 9 and 17 with parsed gene annotation

## macs2

+ bdg - Folder holding the treatment fpkm profiles generated via macs2 for A549 samples and two simulated replicates generted in Roth et al. 

### CHiPR bed files 

+ ATAC_TRIAL1_chipr_peaks_all.bed - CHiPR IDR peaks from first three A549 ATAC-seq replicates on fresh cells

+ Cryo_TRIAL1_chipr_peaks_all.bed - CHiPR IRD peaks from two, cryo-preserverd, replicate A549, ATAC-seq expeirments 

+ Cryo_TRIAL2_chipr_peaks_all.bed - CHiPR IRD peaks from three, cryo-preserverd, replicate A549, ATAC-seq expeirments 

## MISC

Miscellaneous data generated during and for analysis

### SYNTHETIC 

+ Fpkm counts from various synthetic replicates 

## SRA

SRA input genrated for submission to NCBI's sequence read archive. 