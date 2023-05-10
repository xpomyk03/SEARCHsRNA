# It is required to install this packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")

install.packages("seqinr")

install.packages("dplyr")