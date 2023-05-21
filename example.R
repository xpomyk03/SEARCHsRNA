# Example script for running SEARCHsRNA tool
# Uploading of packages
library(Rsamtools)
library(seqinr)
library(dplyr)

# Loading all necessary functions for SEARCHsRNA tool
source('~/SEARCHsRNA/functions.R')

# Setting the path to the folder containing the data - BAM, FASTA, GFF3
path_of_files <- '~/path_to_my_data' 


#### Parameters of function search_sRNA:
# type_of_data:                           'Stranded' for stranded BAM data/ 'ReverslyStranded' for Reversly stranded BAM data [default: 'Stranded']
# threshold_coverage_sRNA                  Set this value [] as the minimum coverage requirement of the resulting sRNAs [default: 0]
# min_length_of_sRNA                       Set this value [pb] as the minimum length requirement of the resulting sRNAs [default: 40]
# threshold_coverage_steepness             Set this value [number of reads] to set the value to detect changes in coverage and find transcripts (a higher value searches only for transcripts with high steepness coverage between  the 5'/3' ends and introns) [default: NULL]
# threshold_coverage_min                   Set this value [number of minimal alignment reads] to set the value of minimum coverage for 3' and 5' ends (a higher value searches only the 5' and 3' ends with higher coverage) [default: NULL]
# threshold_gap_transcripts                Set this value [pb] to determine the minimum gap that is allowed between 2 transcripts (if the gap is smaller, the transcripts will be joined together) [default: NULL]


#Example for manual settings of parametters
search_sRNA(path_of_files, type_of_data = 'Stranded', threshold_coverage_sRNA = 10, min_length_of_sRNA = 50, threshold_coverage_steepness = 5, threshold_coverage_min = 2, threshold_gap_transcripts = 50)

#Example for automatic settings of parametters
search_sRNA(path_of_files)
