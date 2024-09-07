aggregate_counts <- function(cntTable_complete.list) {
  
  # cntTable_complete.list should be a list of counts files
  
  
  
  # Read 1 count table and keep only the gene column. In the FOR loop below, extend by adding counts from each sample. 
  count_data <- read.csv(cntTable_complete.list[1], header=TRUE, sep="")
  count_data <- count_data[1]
  
  
  for (count_table_file in cntTable_complete.list) {
    
    # Load individual table temporarily. 
    count_data_temp <- read.csv(count_table_file, header=TRUE, sep="")
    
    # Check whether the GENE/TE names column matches that of the main count_data table.  
    column_1_match <- all.equal(count_data_temp[1], count_data[1])
    
    if (column_1_match == TRUE) {
      
      # Add the temporary counts data to the final counts object.
      count_data <- data.frame(count_data, count_data_temp[2])
      
      # Remove count_data_temp since it is not needed
      rm(count_data_temp)
      
    } else {
      
      count_data <- c('Count table gene names do not match in some samples')
      break
      
    } # End if else.
    
  } #End FOR loop
  
  
  
  # Output the final count data table.
  return(count_data)
  
  
  
  
} # END FUNCTION

cleanup_and_filter_counts <- function(count_data, filter, min_counts_per_sample, fraction_of_samples) {
  
  # FUNCTION INFO:
  # CLEANUP: REMOVE GENE/TRANSCRIPT VERSION AND COMBINE COUNTS FROM THE SAME GENE
  # FILTER: LOWLY EXPRESSED GENES
  
  
  
  # Assign name to the column with gene/TE names. Name needs to specifically be 'gene.TE' since that is referenced downstream.
  names(count_data)[1] <- 'gene.TE'    

  
  # Remove transcript or ENSEMBL version info. 
  count_data[[1]] <- sub("\\..*", "", count_data[[1]], fixed=FALSE)
  
  
  # Collect same transcripts into one row. NOTE, pseudoautosomal PAR_Y genes will be combined with genes in homologous regions.
  gene_counts_summed <- aggregate(count_data[,-c(1)],
                                  by = list(count_data$gene.TE),
                                  FUN = 'sum')
  
  # Update the counts matrix by moving the gene/TE column to rownames and deleting the original column
  rownames(gene_counts_summed) <- gene_counts_summed[, 1]
  gene_counts_summed <- gene_counts_summed[, c(-1)]
    
    
  # Remove lowly expressed genes, if desired. Filter by CPM (about 10 reads across all samples) to control for differences in library depth.
  if (filter == TRUE) {
    
    # Calculate the number of samples in the counts table
    num_of_samples <- ncol(count_data) - 1
    
    # Define the minimum number of samples that have to meet the CPM threshold
    min.samples <- round(fraction_of_samples * num_of_samples)
    
    # Calculate library sizes for each sample 
    library_sizes <- colSums(gene_counts_summed)
    
    # Calculate the median library size
    MedianLibSize <- median(library_sizes)
    
    # calculate the CPM cutoff for the input libraries (changes with the median library size)
    CPM.Cutoff <- min_counts_per_sample / MedianLibSize*1e6
    
    # Convert counts to CPMs
    CPM_table <- edgeR::cpm(gene_counts_summed, lib.size = library_sizes)
    
    # Make a filter function, where "min.samples" have to exceed "CPM.Cutoff"
    cpm.filter.func <- genefilter::kOverA(min.samples, CPM.Cutoff)
    
    # Bind the filtering function to filterfun
    flist <- genefilter::filterfun(cpm.filter.func)
    
    # Identify genes that pass the filtering threshold
    kept_genes <- genefilter::genefilter(CPM_table, flist)
    
    # Get counts for the genes that pass filtering
    gene_counts_summed_filtered <- gene_counts_summed[kept_genes, ]
    
    # Return the cleaned and filtered counts
    return(gene_counts_summed_filtered)
    
  } else {
    
    # Return the cleaned counts
    return(gene_counts_summed)
    
  }
  
  

  
} # END FUNCTION
