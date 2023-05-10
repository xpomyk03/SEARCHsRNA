# Uploading of packages
library(Rsamtools)
library(seqinr)
library(dplyr)

# Loading all necessary scripts
source('~/DP_sRNA/function.R')

# Set the working directory
setwd('~/DP_sRNA')

# Setting of parameters
type_of_data <- 'ReverseStranded'  # 'Stranded' for stranded BAM data/ 'ReverseStranded' for reverse stranded BAM data
threshold_coverage_transcripts_user <- 10 #Set this value as the minimum coverage requirement of the resulting sRNAs
min_length_of_sRNA_user <- 50 #Set this value as the minimum length requirement of the resulting sRNAs
  # Pro
threshold_coverage_jump_user  <- NULL #Set this value [number of reads] to set the value to detect changes in coverage and find transcripts (a higher value searches only for transcripts with high jump coverage between  the 5'/3' ends and introns)
threshold_coverage_min_user <- NULL #Set this value [number of minimal alignment reads] to set the value of minimum coverage for 3' and 5' ends (a higher value searches only the 5' and 3' ends with higher coverage)
threshold_gap_transkripts_user <- 250 #Set this value [pb] to determine the minimum gap that is allowed between 2 transcripts (if the gap is smaller, the transcripts will be joined together)


# Loading of reference - fasta
fastafile <- dir('.', 'fasta$'); fasta <- read.fasta(fastafile); rm(fastafile);

#Loading of annotation - gff3
gfffile <- dir('.', 'gff3$'); gff <- read.table(gfffile, sep='\t', quote=''); rm(gfffile);

## Look for all files in the current directory that end with 'bam'
bamfiles <- dir('.', 'bam$')

#Function for search sRNA
search_sRNA(bamfiles, gff, fasta, threshold_coverage_transcripts_user, min_length_of_sRNA_user, type_of_data, threshold_coverage_jump_user, threshold_coverage_min_user, threshold_gap_transkripts_user)

