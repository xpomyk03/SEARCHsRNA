# Uploading of packages
library(Rsamtools)
library(seqinr)
library(dplyr)

# Loading all necessary functions for SEARCHsRNA tool
source('~/DP_sRNA_searching/functions.R')

# Setting the path to the folder containing the data
path_of_files <- '~/DP_sRNA_searching' 

#Function for search sRNA
search_sRNA(path_of_files, threshold_coverage_sRNA_user = 10)

