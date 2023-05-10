# FUNKCE Z https://gist.github.com/SamBuckberry/9914246
readBAM <- function(bamFile){
  
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


search_sRNA <- function(bamfiles, gff, fasta, threshold_coverage_transcripts_user, min_length_of_sRNA_user, type_of_data, threshold_coverage_jump_user, threshold_coverage_min_user, threshold_gap_transkripts_user){
  length_of_genom <- length(fasta[[1]])
  
  # Annotation of positive string
  anoted_genes_start_positive <- c()
  anoted_genes_end_positive <- c()
  anoted_genes_start_positive <- gff[,4]
  anoted_genes_end_positive <- gff[,5]
  anoted_genes_start_positive <- anoted_genes_start_positive[(gff[,3] == 'CDS') & (gff[,7] == '+')]
  anoted_genes_end_positive <- anoted_genes_end_positive[(gff[,3] == 'CDS') & (gff[,7] == '+')]
  
  
  #jen pro vizualizaci signal - SMAZAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  anoted_signal_positive <- rep(0, length_of_genom)
  for(k in 1:length(anoted_genes_end_positive)){
    anoted_signal_positive[anoted_genes_start_positive[k]:anoted_genes_end_positive[k]] <- 1
  }
  
  # Annotation of negative string
  anoted_genes_start_negative <- c()
  anoted_genes_end_negative <- c()
  anoted_genes_start_negative <- gff[,4]
  anoted_genes_end_negative <- gff[,5]
  anoted_genes_start_negative <- anoted_genes_start_negative[(gff[,3] == 'CDS') & (gff[,7] == '-')]
  anoted_genes_end_negative <- anoted_genes_end_negative[(gff[,3] == 'CDS') & (gff[,7] == '-')]
  
  #jen pro vizualizaci signal - SMAZAT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  anoted_signal_negative <- rep(0, length_of_genom)
  for(k in 1:length(anoted_genes_end_negative)){
    anoted_signal_negative[anoted_genes_start_negative[k]:anoted_genes_end_negative[k]] <- 1
  }
  rm(k)
  
  
  for (i in 1:length(bamfiles)){
    filename <- bamfiles[i]
    
    ## Load the BAM file
    BAM <- readBAM(filename)
    
    ### Search for potential transcripts
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
    
    rm(list = c('pos', 'strand', 'qwidth'))
    
    # Signal of coverage for positive string
    signal_positive <- rep(0, length_of_genom)
    for (i in 1:length(pos_positive)){
      signal_positive[(pos_positive[i]):(pos_positive[i]+qwidth_positive[i]-1)] <- signal_positive[(pos_positive[i]):(pos_positive[i]+qwidth_positive[i]-1)] + rep(1, qwidth_positive[i])
    }
    signal_positive <- signal_positive[1:length_of_genom]
    
    # Signal of coverage for negative string
    signal_negative <- rep(0, length_of_genom)
    for (i in 1:length(pos_negative)){
      signal_negative[(pos_negative[i]):(pos_negative[i]+qwidth_negative[i]-1)] <- signal_negative[(pos_negative[i]):(pos_negative[i]+qwidth_negative[i]-1)] + rep(1, qwidth_negative[i])
    }
    signal_negative <- signal_negative[1:length_of_genom]
    
    # Coverage of reads
    coverage_signal <- (sum(signal_positive)+sum(signal_negative))/(2*length_of_genom)
    
    # Setting threshold of 'jump' coverage
    if(is.null(threshold_coverage_jump_user)){
      threshold_coverage_jump <- (((-2/700)*coverage_signal) + (195/700))*coverage_signal
    }else{
      threshold_coverage_jump <- threshold_coverage_jump_user
    }
    
    # Setting threshold of minimal coverage
    if(is.null(threshold_coverage_min_user)){
      threshold_coverage_min <- ((-(1/600)*coverage_signal) + (13/60))*coverage_signal
    }else{
      threshold_coverage_min <- threshold_coverage_min_user
    }
    
    # Setting threshold of gap between transcripts
    if(is.null(threshold_gap_transkripts_user)){
      threshold_gap_transkripts <- ((-10/7)*coverage_signal)+151
    }else{
      threshold_gap_transkripts <- threshold_gap_transkripts_user
    }
    
    if(is.null(threshold_coverage_transcripts_user)){
      threshold_coverage_transcripts <- Inf
    }else{
      threshold_coverage_transcripts <- threshold_coverage_transcripts_user
    }
    
    
    
    # POSITIVE STRING
    # Searching for 5' and 3' potential ends of transcripts
    # 5' ends
    position_start_positive <- rep(0, length_of_genom)
    threshold_high <- Inf
    threshold_low <- threshold_coverage_min
    for (k in 1:(length_of_genom-29)){
      if((signal_positive[k+29]-signal_positive[k]) > threshold_coverage_jump){
        if((signal_positive[k+29] < threshold_high) & (signal_positive[k] > threshold_low)){
          a <- signal_positive[k:(k+29)]
          a[a == 0] <- Inf
          position <- which(a == (min(a)))
          position <- position[1]
          position_start_positive[k+position-1] <- 1
          threshold_high <- threshold_coverage_min
          threshold_low <- threshold_coverage_min
        }
      }
      if(signal_positive[k]  < threshold_low){
        threshold_low <- threshold_coverage_min
        threshold_high <- Inf
      }
    }
    
    # 3' ends
    position_end_positive <- rep(0, length_of_genom)
    threshold_high <- Inf
    threshold_low <- threshold_coverage_min
    for (k in (length_of_genom:30)){
      if((signal_positive[k-29]-signal_positive[k]) > threshold_coverage_jump){
        if((signal_positive[k-29] < threshold_high) & (signal_positive[k] > threshold_low)){
          a <- signal_positive[(k-29):k]
          a[a == 0] <- Inf
          position <- which(a == (min(a)))
          position <- position[length(position)]
          position_end_positive[k-30+position] <- 1
          threshold_high <- threshold_coverage_min
          threshold_low <- threshold_coverage_min
        }
      }
      if(signal_positive[k]  < threshold_low){
        threshold_low <- threshold_coverage_min
        threshold_high <- Inf
      }
    }
    
    rm(list = c('threshold_high', 'threshold_low', 'position', 'k', 'a'))
    
    
    
    range <-  c(20000:35000) #10000-15000 40000-50000
    plot(signal_positive[range], type='l', lwd=2, col='red')  #22400:23300
    par(new = TRUE)
    plot(signal_negative[range], type='l', lwd=2, col='black')  #22400:23300
    par(new = TRUE)
    plot(position_start_positive[range], type='l', lwd=2, col='blue')
    par(new = TRUE)
    plot(position_end_positive[range], type='l', lwd=2, col='green')
    par(new = TRUE)
    plot(anoted_signal_positive[range], type='l', lwd=2, col='orange')
    par(new = TRUE)
    plot(anoted_signal_negative[range], type='l', lwd=2, col='brown')
    
    
    
    # Deleting of false ends
    position_start_end_positive <- position_start_positive
    position_start_end_positive[position_end_positive == 1] <- 2
    pamet <- 0
    for (k in 1:(length_of_genom)){
      if ((pamet == 1) & (position_start_end_positive[k] == 1)){
        index_of_start_2 <- k
        s <- 0
        e <- k+1
        while(s == 0){
          if(position_start_end_positive[e] == 2){
            index_of_end <- e
            s <- 2
          }else{
            e <- e+1
          }
        }
        is_null_1 <- min(signal_positive[index_of_start_1:index_of_end])
        if(is_null_1 < 0.5*threshold_coverage_min){
          position_start_end_positive[index_of_start_1] <- 0
          index_of_start_1 <- index_of_start_2
        }else{
          position_start_end_positive[index_of_start_2] <- 0
        }
      }else if((pamet == 0) & (position_start_end_positive[k] == 1)){
        p <- 0
        p2 <- k+1
        while((p == 0) & (p2 < (length_of_genom+1))){
          if(position_start_end_positive[p2] == 2){
            pamet <- 1
            index_of_start_1 <- k
            p <- 1
          }else if(p2 == length_of_genom){
            position_start_end_positive[k] <- 0
            p2 <- p2+1
          }else{
            p2 <- p2+1
          }
        }
        if(k == length_of_genom){
          position_start_end_positive[k] <- 0
        }
      }else if((pamet == 1) & (position_start_end_positive[k] == 2)){
        pamet <- 0
      }
    }
    
    pamet <- 0
    for (k in (length_of_genom:1)){
      if ((pamet == 2) & (position_start_end_positive[k] == 2)){
        index_of_end_2 <- k
        s <- 0
        e <- k-1
        while(s == 0){
          if(position_start_end_positive[e] == 1){
            index_of_start <- e
            s <- 2
          }else{
            e <- e-1
          }
        }
        is_null_1 <- min(signal_positive[index_of_start:index_of_end_1])
        if(is_null_1 < 0.5*threshold_coverage_min){
          position_start_end_positive[index_of_end_1] <- 0
          index_of_end_1 <- index_of_end_2
        }else{
          position_start_end_positive[index_of_end_2] <- 0
        }
      }else if((pamet == 0) & (position_start_end_positive[k] == 2)){
        p <- 0
        p2 <- k-1
        while((p == 0) & (p2 > 0)){
          if(position_start_end_positive[p2] == 1){
            pamet <- 2
            index_of_end_1 <- k
            p <- 1
          }else if(p2 == 1){
            position_start_end_positive[k] <- 0
            p2 <- p2-1
          }else{
            p2 <- p2-1
          }
        }
        if(k == 1){
          position_start_end_positive[k] <- 0
        }
      }else if((pamet == 2) & (position_start_end_positive[k] == 1)){
        pamet <- 0
      }
    }
    
    position_start_positive <- position_start_positive * position_start_end_positive
    position_end_positive <- (position_end_positive * position_start_end_positive)/2
    
    rm(list = c('p','p2','pamet', 'index_of_end_1', 'index_of_end_2', 'e', 's', 'index_of_start', 'index_of_start_1', 'index_of_start_2', 'index_of_end', 'is_null_1', 'position_start_end_positive'))
    
    # Extension of 3' and 5' ends in both directions
    start_index_positive <- c()
    end_index_positive <- c()
    
    for(k in 1:length_of_genom){
      if(position_start_positive[k] == 1){
        start_index_positive <- c(start_index_positive, k)
      }else if(position_end_positive[k] == 1){
        end_index_positive <- c(end_index_positive, k)
      }
    }
    
    for(k in 1:(length(start_index_positive))){
      s <- 0
      p1 <- start_index_positive[k]
      while ((s < 1) & (p1 > 1)){
        if (signal_positive[p1-1] > (threshold_coverage_min)){
          start_index_positive[k] <- p1-1
          p1 <- p1 - 1
        }else{
          s <- 1
        }
      }
      e <- 0
      p2 <- end_index_positive[k]
      while ((e < 1) & (p2 < (length_of_genom-1))){
        if (signal_positive[p2+1] > (threshold_coverage_min)){
          end_index_positive[k] <- p2+1
          p2 <- p2 + 1
        }else{
          e <- 1
        }
      }
    }
    
    position_start_positive <- rep(0, length_of_genom)
    position_end_positive <- rep(0, length_of_genom)
    position_start_positive[start_index_positive] <- 1
    position_end_positive[end_index_positive] <- 1
    
    rm(list = c('s', 'k', 'e', 'p1', 'p2'))
    
    
    # Coverage for potential transcripts
    mean_coverage_transcripts_positive <- rep(0, length(start_index_positive))
    for (i in 1:length(start_index_positive)){
      mean_coverage_transcripts_positive[i] <- mean(signal_positive[start_index_positive[i]:end_index_positive[i]])
    }
    
    # Connection of nearby transcripts
    end_index_positive_clear <- end_index_positive
    start_index_positive_clear <- start_index_positive
    for(k in 1:(length(start_index_positive)-1)){
      if(((start_index_positive[k+1] - end_index_positive[k]) < 20)){
        end_index_positive_clear[k] <- 0
        start_index_positive_clear[k+1] <-0
      }
      else if((abs(mean_coverage_transcripts_positive[k+1]-mean_coverage_transcripts_positive[k]) > 0.5*coverage_signal)){
        next
      }else if(((start_index_positive[k+1] - end_index_positive[k]) < threshold_gap_transkripts)){
        end_index_positive_clear[k] <- 0
        start_index_positive_clear[k+1] <-0
      }
    }
    start_index_positive <- start_index_positive[start_index_positive_clear != 0]
    end_index_positive <- end_index_positive[end_index_positive_clear != 0]
    rm(list = c('start_index_positive_clear', 'end_index_positive_clear', 'k', 'i'))
    
    # signal detekovanych transcripts - respektice oznaceni, kde jsou predikovane, jen pro vizualizaci, pak SMAZAT!
    transkript_signal_positive <- rep(0, length_of_genom)
    for (k in 1:length(start_index_positive)){
      transkript_signal_positive[start_index_positive[k]:end_index_positive[k]] <- 0.5
    }
    
    position_start_positive <- rep(0, length_of_genom)
    position_end_positive <- rep(0, length_of_genom)
    position_start_positive[start_index_positive] <- 1
    position_end_positive[end_index_positive] <- 1
    
    
    range <- c(22000:35000)
    plot(signal_positive[range], type='l', lwd=2, col='red')  #22400:23300
    par(new = TRUE)
    plot(signal_negative[range], type='l', lwd=2, col='black')  #22400:23300
    par(new = TRUE)
    plot(position_start_positive[range], type='l', lwd=2, col='blue')
    par(new = TRUE)
    plot(position_end_positive[range], type='l', lwd=2, col='green')
    par(new = TRUE)
    plot(anoted_signal_positive[range], type='l', lwd=2, col='orange')
    par(new = TRUE)
    plot(anoted_signal_negative[range], type='l', lwd=2, col='brown')
    par(new = TRUE)
    plot(transkript_signal_positive[range], type='l', lwd=2, col='pink')
    
    
    # Deleting of potential transcripts by length
    length_of_transcripts_positive <- end_index_positive - start_index_positive + 1
    start_index_positive <- start_index_positive[length_of_transcripts_positive > min_length_of_sRNA_user]
    end_index_positive <- end_index_positive[length_of_transcripts_positive > min_length_of_sRNA_user]
    length_of_transcripts_positive <- length_of_transcripts_positive[length_of_transcripts_positive > min_length_of_sRNA_user]
    
    # Deleting of potential transcripts by annotation
    start_index_positive_without <- start_index_positive
    end_index_positive_without <- end_index_positive
    for(k in 1:length(anoted_genes_end_positive)){
      start_index_positive_without[((start_index_positive < (anoted_genes_end_positive[k]+50)) & (start_index_positive > (anoted_genes_start_positive[k]-50))) | ((end_index_positive < (anoted_genes_end_positive[k]+50)) & (end_index_positive > (anoted_genes_start_positive[k]-50))) | ((start_index_positive < anoted_genes_start_positive[k]) & (end_index_positive > anoted_genes_end_positive[k]))] <- 0
      end_index_positive_without[((start_index_positive < (anoted_genes_end_positive[k]+50)) & (start_index_positive > (anoted_genes_start_positive[k]-50))) | ((end_index_positive < (anoted_genes_end_positive[k]+50)) & (end_index_positive > (anoted_genes_start_positive[k]-50))) | ((start_index_positive < anoted_genes_start_positive[k]) & (end_index_positive > anoted_genes_end_positive[k]))] <- 0
    }
    start_index_positive <- start_index_positive[start_index_positive_without != 0]
    end_index_positive <- end_index_positive[end_index_positive_without != 0]
    
    rm(list = c('start_index_positive_without', 'end_index_positive_without'))
    
    # Deleting of potential transcripts by mean coverage of transcript
    mean_coverage_transcripts_positive <- rep(0, length(start_index_positive))
    for (i in 1:length(start_index_positive)){
      mean_coverage_transcripts_positive[i] <- mean(signal_positive[start_index_positive[i]:end_index_positive[i]])
    }
    start_index_positive[mean_coverage_transcripts_positive < threshold_coverage_transcripts] <- 0
    end_index_positive[mean_coverage_transcripts_positive < threshold_coverage_transcripts] <- 0
    mean_coverage_transcripts_positive[mean_coverage_transcripts_positive < threshold_coverage_transcripts] <- 0
    start_index_positive <- start_index_positive[start_index_positive != 0]
    end_index_positive <- end_index_positive[end_index_positive != 0]
    mean_coverage_transcripts_positive <- mean_coverage_transcripts_positive[mean_coverage_transcripts_positive != 0]
    
    
    # signal detekovanych transcripts - respektice oznaceni, kde jsou predikovane ---- pro zobrazeni?
    transkript_signal_positive <- rep(0, length_of_genom)
    for (k in 1:length(start_index_positive)){
      transkript_signal_positive[start_index_positive[k]:end_index_positive[k]] <- 1
    }
    position_start_positive <- rep(0, length_of_genom)
    position_end_positive <- rep(0, length_of_genom)
    position_start_positive[start_index_positive] <- 1
    position_end_positive[end_index_positive] <- 1
    
    
    range <-  c(38000:45000) #22000:35000
    plot(signal_positive[range], type='l', lwd=2, col='red')  #22400:23300
    par(new = TRUE)
    plot(signal_negative[range], type='l', lwd=2, col='black')  #22400:23300
    par(new = TRUE)
    plot(anoted_signal_positive[range], type='l', lwd=2, col='orange')
    par(new = TRUE)
    plot(anoted_signal_negative[range], type='l', lwd=2, col='brown')
    par(new = TRUE)
    plot(transkript_signal_positive[range], type='l', lwd=2, col='pink')
    
    
    
    # NEGATIVE STRING
    # Searching for 5' and 3' potential ends of transcripts
    # 5' ends
    position_start_negative <- rep(0, length_of_genom)
    threshold_high <- Inf
    threshold_low <- threshold_coverage_min
    for (k in 1:(length_of_genom-29)){
      if((signal_negative[k+29]-signal_negative[k]) > threshold_coverage_jump){
        if((signal_negative[k+29] < threshold_high) & (signal_negative[k] > threshold_low)){
          a <- signal_negative[k:(k+29)]
          a[a == 0] <- Inf
          position <- which(a == (min(a)))
          position <- position[1]
          position_start_negative[k+position-1] <- 1
          threshold_high <- threshold_coverage_min
          threshold_low <- threshold_coverage_min
        }
      }
      if(signal_negative[k]  < threshold_low){
        threshold_low <- threshold_coverage_min
        threshold_high <- Inf
      }
    }
    
    # 3' ends
    position_end_negative <- rep(0, length_of_genom)
    threshold_high <- Inf
    threshold_low <- threshold_coverage_min
    for (k in (length_of_genom:30)){
      if((signal_negative[k-29]-signal_negative[k]) > threshold_coverage_jump){
        if((signal_negative[k-29] < threshold_high) & (signal_negative[k] > threshold_low)){
          a <- signal_negative[(k-29):k]
          a[a == 0] <- Inf
          position <- which(a == (min(a)))
          position <- position[length(position)]
          position_end_negative[k-30+position] <- 1
          threshold_high <- threshold_coverage_min
          threshold_low <- threshold_coverage_min
        }
      }
      if(signal_negative[k]  < threshold_low){
        threshold_low <- threshold_coverage_min
        threshold_high <- Inf
      }
    }
    
    rm(list = c('threshold_high', 'threshold_low', 'position', 'k', 'a'))
    
    # Deleting of false ends
    position_start_end_negative <- position_start_negative
    position_start_end_negative[position_end_negative == 1] <- 2
    pamet <- 0
    for (k in 1:(length_of_genom)){
      if ((pamet == 1) & (position_start_end_negative[k] == 1)){
        index_of_start_2 <- k
        s <- 0
        e <- k+1
        while(s == 0){
          if(position_start_end_negative[e] == 2){
            index_of_end <- e
            s <- 2
          }else{
            e <- e+1
          }
        }
        is_null_1 <- min(signal_negative[index_of_start_1:index_of_end])
        if(is_null_1 < 0.5*threshold_coverage_min){
          position_start_end_negative[index_of_start_1] <- 0
          index_of_start_1 <- index_of_start_2
        }else{
          position_start_end_negative[index_of_start_2] <- 0
        }
      }else if((pamet == 0) & (position_start_end_negative[k] == 1)){
        p <- 0
        p2 <- k+1
        while((p == 0) & (p2 < (length_of_genom+1))){
          if(position_start_end_negative[p2] == 2){
            pamet <- 1
            index_of_start_1 <- k
            p <- 1
          }else if(p2 == length_of_genom){
            position_start_end_negative[k] <- 0
            p2 <- p2+1
          }else{
            p2 <- p2+1
          }
        }
        if(k == length_of_genom){
          position_start_end_positive[k] <- 0
        }
      }else if((pamet == 1) & (position_start_end_negative[k] == 2)){
        pamet <- 0
      }
    }
    
    pamet <- 0
    for (k in (length_of_genom:1)){
      if ((pamet == 2) & (position_start_end_negative[k] == 2)){
        index_of_end_2 <- k
        s <- 0
        e <- k-1
        while(s == 0){
          if(position_start_end_negative[e] == 1){
            index_of_start <- e
            s <- 2
          }else{
            e <- e-1
          }
        }
        is_null_1 <- min(signal_negative[index_of_start:index_of_end_1])
        if(is_null_1 < 0.5*threshold_coverage_min){
          position_start_end_negative[index_of_end_1] <- 0
          index_of_end_1 <- index_of_end_2
        }else{
          position_start_end_negative[index_of_end_2] <- 0
        }
      }else if((pamet == 0) & (position_start_end_negative[k] == 2)){
        p <- 0
        p2 <- k-1
        while((p == 0) & (p2 > 0)){
          if(position_start_end_negative[p2] == 1){
            pamet <- 2
            index_of_end_1 <- k
            p <- 1
          }else if(p2 == 1){
            position_start_end_negative[k] <- 0
            p2 <- p2-1
          }else{
            p2 <- p2-1
          }
        }
        if(k == 1){
          position_start_end_positive[k] <- 0
        }
      }else if((pamet == 2) & (position_start_end_negative[k] == 1)){
        pamet <- 0
      }
    }
    
    position_start_negative <- position_start_negative * position_start_end_negative
    position_end_negative <- (position_end_negative * position_start_end_negative)/2
    
    rm(list = c('p','p2','pamet', 'index_of_end_1', 'index_of_end_2', 'e', 's', 'index_of_start', 'index_of_start_1', 'index_of_start_2', 'index_of_end', 'is_null_1', 'position_start_end_negative'))
    
    # Extension of 3' and 5' ends in both directions
    start_index_negative <- c()
    end_index_negative <- c()
    for(k in 1:length_of_genom){
      if(position_start_negative[k] == 1){
        start_index_negative <- c(start_index_negative, k)
      }else if(position_end_negative[k] == 1){
        end_index_negative <- c(end_index_negative, k)
      }
    }
    
    for(k in 1:(length(start_index_negative))){
      s <- 0
      p1 <- start_index_negative[k]
      while ((s < 1) & (p1 > 1)){
        if (signal_negative[p1-1] > (threshold_coverage_min)){
          start_index_negative[k] <- p1-1
          p1 <- p1 - 1
        }else{
          s <- 1
        }
      }
      e <- 0
      p2 <- end_index_negative[k]
      while ((e < 1) & (p2 < (length_of_genom-1))){
        if (signal_negative[p2+1] > (threshold_coverage_min)){
          end_index_negative[k] <- p2+1
          p2 <- p2 + 1
        }else{
          e <- 1
        }
      }
    }
    
    position_start_negative <- rep(0, length_of_genom)
    position_end_negative <- rep(0, length_of_genom)
    position_start_negative[start_index_negative] <- 1
    position_end_negative[end_index_negative] <- 1
    
    rm(list = c('s', 'k', 'e', 'p1', 'p2'))
    
    # Coverage for potential transcripts
    mean_coverage_transcripts_negative <- rep(0, length(start_index_negative))
    for (i in 1:length(start_index_negative)){
      mean_coverage_transcripts_negative[i] <- mean(signal_negative[start_index_negative[i]:end_index_negative[i]])
    }
    
    # Connection of nearby transcripts
    end_index_negative_clear <- end_index_negative
    start_index_negative_clear <- start_index_negative
    for(k in 1:(length(start_index_negative)-1)){
      if(((start_index_negative[k+1] - end_index_negative[k]) < 20)){
        end_index_negative_clear[k] <- 0
        start_index_negative_clear[k+1] <-0
      }
      else if((abs(mean_coverage_transcripts_negative[k+1]-mean_coverage_transcripts_negative[k]) > 0.5*coverage_signal)){
        next
      }else if(((start_index_negative[k+1] - end_index_negative[k]) < threshold_gap_transkripts)){
        end_index_negative_clear[k] <- 0
        start_index_negative_clear[k+1] <-0
      }
    }
    
    start_index_negative <- start_index_negative[start_index_negative_clear != 0]
    end_index_negative <- end_index_negative[end_index_negative_clear != 0]
    
    rm(list = c('start_index_negative_clear', 'end_index_negative_clear', 'k', 'i'))
    
    # signal detekovanych transcripts - respektice oznaceni, kde jsou predikovane, jen pro vizualizaci, pak SMAZAT!
    transkript_signal_negative <- rep(0, length_of_genom)
    for (k in 1:length(start_index_negative)){
      transkript_signal_negative[start_index_negative[k]:end_index_negative[k]] <- 0.5
    }
    
    position_start_negative <- rep(0, length_of_genom)
    position_end_negative <- rep(0, length_of_genom)
    position_start_negative[start_index_negative] <- 1
    position_end_negative[end_index_negative] <- 1
    
    
    
    range <-  c(40000:50000)
    plot(signal_positive[range], type='l', lwd=2, col='red')  #22400:23300
    par(new = TRUE)
    plot(signal_negative[range], type='l', lwd=2, col='black')  #22400:23300
    par(new = TRUE)
    plot(position_start_negative[range], type='l', lwd=2, col='blue')
    par(new = TRUE)
    plot(position_end_negative[range], type='l', lwd=2, col='green')
    par(new = TRUE)
    plot(anoted_signal_positive[range], type='l', lwd=2, col='orange')
    par(new = TRUE)
    plot(anoted_signal_negative[range], type='l', lwd=2, col='brown')
    par(new = TRUE)
    plot(transkript_signal_negative[range], type='l', lwd=2, col='pink')
    
    
    
    # Deleting of potential transcripts by length
    length_of_transcripts_negative <- end_index_negative - start_index_negative + 1
    start_index_negative <- start_index_negative[length_of_transcripts_negative > min_length_of_sRNA_user]
    end_index_negative <- end_index_negative[length_of_transcripts_negative > min_length_of_sRNA_user]
    length_of_transcripts_negative <- length_of_transcripts_negative[length_of_transcripts_negative > min_length_of_sRNA_user]
    
    # Deleting of potential transcripts by annotation
    start_index_negative_without <- start_index_negative
    end_index_negative_without <- end_index_negative
    
    for(k in 1:length(anoted_genes_end_negative)){
      start_index_negative_without[((start_index_negative < (anoted_genes_end_negative[k]+50)) & (start_index_negative > (anoted_genes_start_negative[k]-50))) | ((end_index_negative < (anoted_genes_end_negative[k]+50)) & (end_index_negative > (anoted_genes_start_negative[k]-50))) | ((start_index_negative < anoted_genes_start_negative[k]) & (end_index_negative > anoted_genes_end_negative[k]))] <- 0
      end_index_negative_without[((start_index_negative < (anoted_genes_end_negative[k]+50)) & (start_index_negative > (anoted_genes_start_negative[k]-50))) | ((end_index_negative < (anoted_genes_end_negative[k]+50)) & (end_index_negative > (anoted_genes_start_negative[k]-50))) | ((start_index_negative < anoted_genes_start_negative[k]) & (end_index_negative > anoted_genes_end_negative[k]))] <- 0
    }
    
    start_index_negative <- start_index_negative[start_index_negative_without != 0]
    end_index_negative <- end_index_negative[end_index_negative_without != 0]
    
    rm(list = c('start_index_negative_without', 'end_index_negative_without'))
    
    
    # Deleting of potential transcripts by mean coverage of transcript
    mean_coverage_transcripts_negative <- rep(0, length(start_index_negative))
    for (i in 1:length(start_index_negative)){
      mean_coverage_transcripts_negative[i] <- mean(signal_negative[start_index_negative[i]:end_index_negative[i]])
    }
    start_index_negative[mean_coverage_transcripts_negative < threshold_coverage_transcripts] <- 0
    end_index_negative[mean_coverage_transcripts_negative < threshold_coverage_transcripts] <- 0
    mean_coverage_transcripts_negative[mean_coverage_transcripts_negative < threshold_coverage_transcripts] <- 0
    start_index_negative <- start_index_negative[start_index_negative != 0]
    end_index_negative <- end_index_negative[end_index_negative != 0]
    mean_coverage_transcripts_negative <- mean_coverage_transcripts_negative[mean_coverage_transcripts_negative != 0]
    
    # signal detekovanych transcripts - respektice oznaceni, kde jsou predikovane ---- pro zobrazeni?
    transkript_signal_negative <- rep(0, length_of_genom)
    for (k in 1:length(start_index_negative)){
      transkript_signal_negative[start_index_negative[k]:end_index_negative[k]] <- 1
    }
    
    position_start_negative <- rep(0, length_of_genom)
    position_end_negative <- rep(0, length_of_genom)
    position_start_negative[start_index_negative] <- 1
    position_end_negative[end_index_negative] <- 1
    
    
    range <-  c(50000:60000) #10000-15000 40000-50000
    plot(signal_positive[range], type='l', lwd=2, col='red')  #22400:23300
    par(new = TRUE)
    plot(signal_negative[range], type='l', lwd=2, col='black')  #22400:23300
    par(new = TRUE)
    plot(anoted_signal_positive[range], type='l', lwd=2, col='orange')
    par(new = TRUE)
    plot(anoted_signal_negative[range], type='l', lwd=2, col='brown')
    par(new = TRUE)
    plot(transkript_signal_negative[range], type='l', lwd=2, col='pink')
    par(new = TRUE)
    plot(transkript_signal_positive[range], type='l', lwd=2, col='blue')
    
    
    # Exporting the results to CSV 
    data <- data.frame(matrix(nrow = 0, ncol = 7))
    columns = c('start', 'stop', 'length', 'strand', 'mean_coverage', 'type', 'gene')
    colnames(data) = columns
    
    length_of_transcripts_positive <- end_index_positive - start_index_positive + 1
    length_of_transcripts_negative <- end_index_negative - start_index_negative + 1
    
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
    data$mean_coverage <- round(as.numeric(data$mean_coverage))
    
    data <- data %>% 
      as.data.frame() %>% 
      arrange(start)
    
    name <-  strsplit(filename, split = "")
    name <- (name[[1]])
    name <- name[1:(length(name)-4)]
    name <- paste(name, collapse='')
    name2 <- '_sRNA.csv'
    name3 <- paste(name, name2, sep = '')
    name2 <- '_positive_sRNA.txt'
    name4 <- paste(name, name2, sep = '')
    name2 <- '_negative_sRNA.txt'
    name5 <- paste(name, name2, sep = '')
    rm(name, name2)
    
    write.table(data, name3, sep = ';', col.names = TRUE, row.names = FALSE) # export
    write.table(transkript_signal_positive, file = name4, sep = '', row.names = FALSE, col.names = FALSE)
    write.table(transkript_signal_negative, file = name5, sep = '', row.names = FALSE, col.names = FALSE)
  }
}

