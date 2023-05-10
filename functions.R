readBAM <- function(bamFile){
  # Function from: https://gist.github.com/SamBuckberry/9914246
  bam <- scanBam(bamFile)
  
  # A function for collapsing the list of lists into a single list
  # as per the Rsamtools vignette
  .unlist <- function (x){
    x1 <- x[[1L]]
    if (is.factor(x1)){
      structure(unlist(x), class = "factor", levels = levels(x1))
    } else {
      do.call(c, x)
    }
  }
  
  bam_field <- names(bam[[1]])
  
  list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
  
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  
  #return a list that can be called as a data frame
  return(bam_df)
}


get_annotation <- function(gff){
  # Annotation of positive string
  genes_start_positive <- c()
  genes_end_positive <- c()
  genes_start_positive <- gff[,4]
  genes_end_positive <- gff[,5]
  genes_start_positive <- genes_start_positive[(gff[,3] == 'CDS') & (gff[,7] == '+')]
  genes_end_positive <- genes_end_positive[(gff[,3] == 'CDS') & (gff[,7] == '+')]
  
  # Annotation of negative string
  genes_start_negative <- c()
  genes_end_negative <- c()
  genes_start_negative <- gff[,4]
  genes_end_negative <- gff[,5]
  genes_start_negative <- genes_start_negative[(gff[,3] == 'CDS') & (gff[,7] == '-')]
  genes_end_negative <- genes_end_negative[(gff[,3] == 'CDS') & (gff[,7] == '-')]
  
  return(list("genes_start_positive" = genes_start_positive, "genes_end_positive" = genes_end_positive, "genes_start_negative" = genes_start_negative, "genes_end_negative" = genes_end_negative))
}


preparing_of_reads <- function(filename, type_of_data){
  # Load the BAM file
  BAM <- readBAM(filename)
  
  # Get information from BAM file
  mapq <- BAM@listData[['mapq']]
  qwidth <- BAM@listData[['qwidth']]
  strand <- BAM@listData[['strand']]
  pos <- BAM@listData[['pos']]
  
  # Removing of low quality reads
  pos <- pos[mapq == 255]
  qwidth <- qwidth[mapq == 255]
  strand <- strand[mapq == 255]
  rm(mapq)
  
  # Removing of really short and long reads
  median_of_length_read <- median(qwidth)
  pos <- pos[(qwidth > 20) & (qwidth < (median_of_length_read+1))]
  strand <- strand[(qwidth > 20) & (qwidth < (median_of_length_read+1))]
  qwidth <- qwidth[(qwidth > 20) & (qwidth < (median_of_length_read+1))]
  
  # Sort the reads for + and - string
  # Type of BAM data setting
  if (type_of_data == 'Stranded'){
    pos_positive <- pos[strand == '+']
    qwidth_positive <- qwidth[strand == '+']
    pos_negative <- pos[strand == '-']
    qwidth_negative <- qwidth[strand == '-']
  }else if(type_of_data == 'ReverseStranded'){
    pos_positive <- pos[strand == '-']
    qwidth_positive <- qwidth[strand == '-']
    pos_negative <- pos[strand == '+']
    qwidth_negative <- qwidth[strand == '+']
  }
  return(list("pos_positive" = pos_positive, "qwidth_positive" = qwidth_positive, "pos_negative" = pos_negative, "qwidth_negative" = qwidth_negative))
}

search_transcripts <- function(signal, length_of_genome, annotation_genes_start, annotation_genes_end, coverage_signal, threshold_coverage_min, threshold_coverage_steepness, threshold_gap_transcripts, min_length_of_sRNA, threshold_coverage_transcripts){
  ## Searching for 5' and 3' potential ends of transcripts
  # 5' ends
  position_start <- rep(0, length_of_genome)
  threshold_high <- Inf
  threshold_low <- threshold_coverage_min
  for (k in 1:(length_of_genome-29)){
    if((signal[k+29]-signal[k]) > threshold_coverage_steepness){
      if((signal[k+29] < threshold_high) & (signal[k] > threshold_low)){
        a <- signal[k:(k+29)]
        a[a == 0] <- Inf
        position <- which(a == (min(a)))
        position <- position[1]
        position_start[k+position-1] <- 1
        threshold_high <- threshold_coverage_min
        threshold_low <- threshold_coverage_min
      }
    }
    if(signal[k]  < threshold_low){
      threshold_low <- threshold_coverage_min
      threshold_high <- Inf
    }
  }
  
  # 3' ends
  position_end <- rep(0, length_of_genome)
  threshold_high <- Inf
  threshold_low <- threshold_coverage_min
  for (k in (length_of_genome:30)){
    if((signal[k-29]-signal[k]) > threshold_coverage_steepness){
      if((signal[k-29] < threshold_high) & (signal[k] > threshold_low)){
        a <- signal[(k-29):k]
        a[a == 0] <- Inf
        position <- which(a == (min(a)))
        position <- position[length(position)]
        position_end[k-30+position] <- 1
        threshold_high <- threshold_coverage_min
        threshold_low <- threshold_coverage_min
      }
    }
    if(signal[k]  < threshold_low){
      threshold_low <- threshold_coverage_min
      threshold_high <- Inf
    }
  }
  
  rm(list = c('threshold_high', 'threshold_low', 'position', 'k', 'a'))
  
  # Deleting of false ends
  position_start_end <- position_start
  position_start_end[position_end == 1] <- 2
  pamet <- 0
  for (k in 1:(length_of_genome)){
    if ((pamet == 1) & (position_start_end[k] == 1)){
      index_of_start_2 <- k
      s <- 0
      e <- k+1
      while(s == 0){
        if(position_start_end[e] == 2){
          index_of_end <- e
          s <- 2
        }else{
          e <- e+1
        }
      }
      is_null_1 <- min(signal[index_of_start_1:index_of_end])
      if(is_null_1 < 0.5*threshold_coverage_min){
        position_start_end[index_of_start_1] <- 0
        index_of_start_1 <- index_of_start_2
      }else{
        position_start_end[index_of_start_2] <- 0
      }
    }else if((pamet == 0) & (position_start_end[k] == 1)){
      p <- 0
      p2 <- k+1
      while((p == 0) & (p2 < (length_of_genome+1))){
        if(position_start_end[p2] == 2){
          pamet <- 1
          index_of_start_1 <- k
          p <- 1
        }else if(p2 == length_of_genome){
          position_start_end[k] <- 0
          p2 <- p2+1
        }else{
          p2 <- p2+1
        }
      }
      if(k == length_of_genome){
        position_start_end[k] <- 0
      }
    }else if((pamet == 1) & (position_start_end[k] == 2)){
      pamet <- 0
    }
  }
  
  pamet <- 0
  for (k in (length_of_genome:1)){
    if ((pamet == 2) & (position_start_end[k] == 2)){
      index_of_end_2 <- k
      s <- 0
      e <- k-1
      while(s == 0){
        if(position_start_end[e] == 1){
          index_of_start <- e
          s <- 2
        }else{
          e <- e-1
        }
      }
      is_null_1 <- min(signal[index_of_start:index_of_end_1])
      if(is_null_1 < 0.5*threshold_coverage_min){
        position_start_end[index_of_end_1] <- 0
        index_of_end_1 <- index_of_end_2
      }else{
        position_start_end[index_of_end_2] <- 0
      }
    }else if((pamet == 0) & (position_start_end[k] == 2)){
      p <- 0
      p2 <- k-1
      while((p == 0) & (p2 > 0)){
        if(position_start_end[p2] == 1){
          pamet <- 2
          index_of_end_1 <- k
          p <- 1
        }else if(p2 == 1){
          position_start_end[k] <- 0
          p2 <- p2-1
        }else{
          p2 <- p2-1
        }
      }
      if(k == 1){
        position_start_end[k] <- 0
      }
    }else if((pamet == 2) & (position_start_end[k] == 1)){
      pamet <- 0
    }
  }
  
  position_start <- position_start * position_start_end
  position_end <- (position_end * position_start_end)/2
  
  rm(list = c('p','p2','pamet', 'index_of_end_1', 'index_of_end_2', 'e', 's', 'index_of_start', 'index_of_start_1', 'index_of_start_2', 'index_of_end', 'is_null_1', 'position_start_end'))
  
  # Extension of 3' and 5' ends in both directions
  start_index <- c()
  end_index <- c()
  
  for(k in 1:length_of_genome){
    if(position_start[k] == 1){
      start_index <- c(start_index, k)
    }else if(position_end[k] == 1){
      end_index <- c(end_index, k)
    }
  }
  
  for(k in 1:(length(start_index))){
    s <- 0
    p1 <- start_index[k]
    while ((s < 1) & (p1 > 1)){
      if (signal[p1-1] > (threshold_coverage_min)){
        start_index[k] <- p1-1
        p1 <- p1 - 1
      }else{
        s <- 1
      }
    }
    e <- 0
    p2 <- end_index[k]
    while ((e < 1) & (p2 < (length_of_genome-1))){
      if (signal[p2+1] > (threshold_coverage_min)){
        end_index[k] <- p2+1
        p2 <- p2 + 1
      }else{
        e <- 1
      }
    }
  }
  
  position_start <- rep(0, length_of_genome)
  position_end <- rep(0, length_of_genome)
  position_start[start_index] <- 1
  position_end[end_index] <- 1
  
  rm(list = c('s', 'k', 'e', 'p1', 'p2'))
  
  
  # Coverage for potential transcripts
  mean_coverage_transcripts <- rep(0, length(start_index))
  for (i in 1:length(start_index)){
    mean_coverage_transcripts[i] <- mean(signal[start_index[i]:end_index[i]])
  }
  
  # Connection of nearby transcripts
  start_index_clear <- start_index
  end_index_clear <- end_index
  for(k in 1:(length(start_index)-1)){
    if(((start_index[k+1] - end_index[k]) < 20)){
      end_index_clear[k] <- 0
      start_index_clear[k+1] <-0
    }
    else if((abs(mean_coverage_transcripts[k+1]-mean_coverage_transcripts[k]) > 0.5*coverage_signal)){
      next
    }else if(((start_index[k+1] - end_index[k]) < threshold_gap_transcripts)){
      end_index_clear[k] <- 0
      start_index_clear[k+1] <-0
    }
  }
  start_index <- start_index[start_index_clear != 0]
  end_index <- end_index[end_index_clear != 0]
  
  rm(list = c('start_index_clear', 'end_index_clear', 'k', 'i'))
  
  
  # Deleting of potential transcripts by minimal length of searched sRNA
  length_of_transcripts <- end_index - start_index + 1
  start_index <- start_index[length_of_transcripts > min_length_of_sRNA]
  end_index <- end_index[length_of_transcripts > min_length_of_sRNA]
  length_of_transcripts <- length_of_transcripts[length_of_transcripts > min_length_of_sRNA]
  
  # Deleting of potential transcripts by annotation
  start_index_without <- start_index
  end_index_without <- end_index
  for(k in 1:length(annotation_genes_end)){
    start_index_without[((start_index < (annotation_genes_end[k]+50)) & (start_index > (annotation_genes_start[k]-50))) | ((end_index < (annotation_genes_end[k]+50)) & (end_index > (annotation_genes_start[k]-50))) | ((start_index < annotation_genes_start[k]) & (end_index > annotation_genes_end[k]))] <- 0
    end_index_without[((start_index < (annotation_genes_end[k]+50)) & (start_index > (annotation_genes_start[k]-50))) | ((end_index < (annotation_genes_end[k]+50)) & (end_index > (annotation_genes_start[k]-50))) | ((start_index < annotation_genes_start[k]) & (end_index > annotation_genes_end[k]))] <- 0
  }
  start_index <- start_index[start_index_without != 0]
  end_index <- end_index[end_index_without != 0]
  rm(list = c('start_index_without', 'end_index_without'))
  
  # Deleting of potential transcripts by mean coverage of transcript
  mean_coverage_transcripts <- rep(0, length(start_index))
  for (i in 1:length(start_index)){
    mean_coverage_transcripts[i] <- mean(signal[start_index[i]:end_index[i]])
  }
  start_index[mean_coverage_transcripts < threshold_coverage_transcripts] <- 0
  end_index[mean_coverage_transcripts < threshold_coverage_transcripts] <- 0
  mean_coverage_transcripts[mean_coverage_transcripts < threshold_coverage_transcripts] <- 0
  start_index <- start_index[start_index != 0]
  end_index <- end_index[end_index != 0]
  mean_coverage_transcripts <- round(mean_coverage_transcripts[mean_coverage_transcripts != 0])
  
  
  # Length of reads
  length_of_transcripts <- end_index - start_index + 1
  
  
  return(list("start_index" = start_index, "end_index" = end_index, "mean_coverage_transcripts" = mean_coverage_transcripts, "length_of_transcripts" = length_of_transcripts))
}



search_sRNA <- function(bamfiles, gff, fasta, threshold_coverage_transcripts_user, min_length_of_sRNA_user, type_of_data, threshold_coverage_steepness_user, threshold_coverage_min_user, threshold_gap_transcripts_user){
  
  # Length of genome from FASTA file
  length_of_genome <- length(fasta[[1]])
  
  # Get positions of annoted genes from annotation FASTA file
  annotation_genes <- get_annotation(gff)
  
  # For-cycle for all BAM files
  for (i in 1:length(bamfiles)){
    filename <- bamfiles[i]
    
    # Get reads from BAM
    reads <- preparing_of_reads(filename, type_of_data) 
    
    # Signal of coverage for positive string
    signal_positive <- rep(0, length_of_genome)
    for (i in 1:length(reads[["pos_positive"]])){
      signal_positive[(reads[["pos_positive"]][i]):(reads[["pos_positive"]][i]+reads[["qwidth_positive"]][i]-1)] <- signal_positive[(reads[["pos_positive"]][i]):(reads[["pos_positive"]][i]+reads[["qwidth_positive"]][i]-1)] + rep(1, reads[["qwidth_positive"]][i])
    }
    signal_positive <- signal_positive[1:length_of_genome]
    
    # Signal of coverage for negative string
    signal_negative <- rep(0, length_of_genome)
    for (i in 1:length(reads[["pos_negative"]])){
      signal_negative[(reads[["pos_negative"]][i]):(reads[["pos_negative"]][i]+reads[["qwidth_negative"]][i]-1)] <- signal_negative[(reads[["pos_negative"]][i]):(reads[["pos_negative"]][i]+reads[["qwidth_negative"]][i]-1)] + rep(1, reads[["qwidth_negative"]][i])
    }
    signal_negative <- signal_negative[1:length_of_genome]
    
    # Coverage of reads
    coverage_signal <- (sum(signal_positive)+sum(signal_negative))/(2*length_of_genome)
    
    ## SETTING THE PARAMETERS
    # Setting threshold of steepness coverage
    if(is.null(threshold_coverage_steepness_user)){
      threshold_coverage_steepness <- (((-2/700)*coverage_signal) + (195/700))*coverage_signal
    }else{
      threshold_coverage_steepness <- threshold_coverage_steepness_user
    }
    
    # Setting threshold of minimal coverage
    if(is.null(threshold_coverage_min_user)){
      threshold_coverage_min <- ((-(1/600)*coverage_signal) + (13/60))*coverage_signal
    }else{
      threshold_coverage_min <- threshold_coverage_min_user
    }
    
    # Setting threshold of gap between transcripts
    if(is.null(threshold_gap_transcripts_user)){
      threshold_gap_transcripts <- ((-10/7)*coverage_signal)+151
    }else{
      threshold_gap_transcripts <- threshold_gap_transcripts_user
    }
    
    # Setting threshold of minimum coverage of final sRNA
    if(is.null(threshold_coverage_transcripts_user)){
      threshold_coverage_transcripts <- Inf
    }else{
      threshold_coverage_transcripts <- threshold_coverage_transcripts_user
    }
    
    # Setting threshold of minimum length of final sRNA
    if(is.null(min_length_of_sRNA_user)){
      min_length_of_sRNA <- 50
    }else{
      min_length_of_sRNA <- min_length_of_sRNA_user
    }
    
    # Searching for sRNA transcripts on positive strand
    positive_transcripts <- search_transcripts(signal_positive, length_of_genome, annotation_genes[["genes_start_positive"]], annotation_genes[["genes_end_positive"]], coverage_signal, threshold_coverage_min, threshold_coverage_steepness, threshold_gap_transcripts, min_length_of_sRNA, threshold_coverage_transcripts)
    
    # Searching for sRNA transcripts on negative strand
    negative_transcripts <- search_transcripts(signal_negative, length_of_genome, annotation_genes[["genes_start_negative"]], annotation_genes[["genes_end_negative"]], coverage_signal, threshold_coverage_min, threshold_coverage_steepness, threshold_gap_transcripts, min_length_of_sRNA, threshold_coverage_transcripts)
    
    # Exporting the results to CSV 
    exporting_CSV(filename, gff, positive_transcripts[["start_index"]], positive_transcripts[["end_index"]], positive_transcripts[["length_of_transcripts"]], positive_transcripts[["mean_coverage_transcripts"]], negative_transcripts[["start_index"]], negative_transcripts[["end_index"]], negative_transcripts[["length_of_transcripts"]], negative_transcripts[["mean_coverage_transcripts"]])
    exporting_signals_TXT(filename, length_of_genome, positive_transcripts, negative_transcripts)
  }
}

exporting_signals_TXT <- function(filename, length_of_genome, positive_transcripts, negative_transcripts){
  
  transkript_signal_positive <- rep(0, length_of_genome)
  for (k in 1:length(positive_transcripts[["start_index"]])){
    transkript_signal_positive[positive_transcripts[["start_index"]][k]:positive_transcripts[["end_index"]][k]] <- 1
  }
  
  transkript_signal_negative <- rep(0, length_of_genome)
  for (k in 1:length(negative_transcripts[["start_index"]])){
    transkript_signal_negative[negative_transcripts[["start_index"]][k]:negative_transcripts[["end_index"]][k]] <- 1
  }
  
  name <- strsplit(filename, split = ""); name <- (name[[1]]); name <- name[1:(length(name)-4)]; name <- paste(name, collapse='')
  nameP <- paste(name, '_positive_sRNA.txt', sep = '')
  nameN <- paste(name, '_negative_sRNA.txt', sep = '')
  
  write.table(transkript_signal_positive, file = nameP, sep = '', row.names = FALSE, col.names = FALSE)
  write.table(transkript_signal_negative, file = nameN, sep = '', row.names = FALSE, col.names = FALSE)
}

exporting_CSV <- function(filename, gff, start_index_positive, end_index_positive, length_of_transcripts_positive, mean_coverage_transcripts_positive, start_index_negative, end_index_negative, length_of_transcripts_negative, mean_coverage_transcripts_negative){
  
  data <- data.frame(matrix(nrow = 0, ncol = 7))
  columns = c('start', 'stop', 'length', 'strand', 'mean_coverage', 'type', 'gene')
  colnames(data) = columns
  
  
  genes <- gff[gff[,3] == 'gene',]
  genes_positive <- genes[genes[,7] == '+',]
  genes_negative <- genes[genes[,7] == '-',]
  
  # positive string
  for(k in 1:length(start_index_positive)){
    list <- genes_negative[((((genes_negative[,4]-1) < start_index_positive[k]) & ((genes_negative[,5]+1) > start_index_positive[k])) | (((genes_negative[,4]-1) < end_index_positive[k]) & ((genes_negative[,5]+1) > end_index_positive[k]))) | ((genes_negative[,4] > start_index_positive[k] & ((genes_negative[,4]) < end_index_positive[k])) | (genes_negative[,5] > start_index_positive[k] & ((genes_negative[,5]) < end_index_positive[k]))),9]
    if(isEmpty(list)){
      type <- 'cis-sRNA'
      gene <- '-'
    }else{
      type <- 'trans-sRNA'
      gene <- list[1]
    }
    new_row <- c(start_index_positive[k], end_index_positive[k], length_of_transcripts_positive[k], '+', mean_coverage_transcripts_positive[k], type, gene)
    data[nrow(data) + 1,] <- new_row  
  }
  
  # negative string
  for(k in 1:length(start_index_negative)){
    list <- genes_positive[((((genes_positive[,4]-1) < start_index_negative[k]) & ((genes_positive[,5]+1) > start_index_negative[k])) | (((genes_positive[,4]-1) < end_index_negative[k]) & ((genes_positive[,5]+1) > end_index_negative[k]))) | ((genes_positive[,4] > start_index_negative[k] & ((genes_positive[,4]) < end_index_negative[k])) | (genes_positive[,5] > start_index_negative[k] & ((genes_positive[,5]) < end_index_negative[k]))),9]
    if(isEmpty(list)){
      type <- 'cis-sRNA'
      gene <- '-'
    }else{
      type <- 'trans-sRNA'
      gene <- list[1]
    }
    new_row <- c(start_index_negative[k], end_index_negative[k], length_of_transcripts_negative[k], '-', mean_coverage_transcripts_negative[k], type, gene)
    data[nrow(data) + 1,] <- new_row  
  }
  
  # sorting the table by start position
  data$start <- strtoi(data$start)
  data$stop <- strtoi(data$stop)
  data$length <- strtoi(data$length)
  data$mean_coverage <- strtoi(data$mean_coverage)
  
  data <- data %>% 
    as.data.frame() %>% 
    arrange(start)
  
  name <- strsplit(filename, split = ""); name <- (name[[1]]); name <- name[1:(length(name)-4)]; name <- paste(name, collapse='')
  name <- paste(name, '_sRNA.csv', sep = '')
  
  write.table(data, name, sep = ';', col.names = TRUE, row.names = FALSE) # export
}