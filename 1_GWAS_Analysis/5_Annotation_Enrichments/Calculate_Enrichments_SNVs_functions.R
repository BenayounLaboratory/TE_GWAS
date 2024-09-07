Annotation_enrichment_sig_vs_background <- function(annotated_background_variants, annotated_significant_variants, annotation_column_name, annotation_label){
  
  

  # Calculate proportions in the significant variant list
     
      # number of snps with the annotation
      sig.in.annotation <- length(which(annotated_significant_variants[, annotation_column_name] == annotation_label))
      
      # number of snps without annotation
      sig.notin.annotation <- nrow(annotated_significant_variants) - sig.in.annotation
      
      # fraction of snps with the annotation
      sig.fraction.annotation <- sig.in.annotation / (sig.in.annotation + sig.notin.annotation)
      
  # Calculate proportions in the background variant list
      
      # number of snps with an annotation
      background.in.annotation <- length(which(all.snvs[, annotation_column_name] == annotation_label))
      
      # number of snps without annotation
      background.notin.annotation <- nrow(all.snvs) - background.in.annotation
      
      # fraction of snps with the annotation
      background.fraction.annotation <- background.in.annotation / (background.in.annotation + background.notin.annotation)
      
  # Collect fractions in a df to use for plot
  both.annotation.fractions <- data.frame(Fraction = c(sig.fraction.annotation, background.fraction.annotation), Group = c('Significant', 'Background'))
      
  # Generate a contingency table as a matrix to run stats
  my.contingency <- as.matrix(data.frame(Outside_Annotation = c(sig.notin.annotation, background.notin.annotation), In_Annotation= c(sig.in.annotation, background.in.annotation)))
  
  # Run Fischer's Exact Test
  my.Fischer.res <- fisher.test(my.contingency)
  
  
  
  # Output 1) the df with fraction of variants with annotation and 2) the Fischer test results
  return(list(both.annotation.fractions, my.Fischer.res))
  
  
} # END FUNCTION

Gene_overlap_enrichment_sig_vs_background <- function(annotated_background_variants, annotated_significant_variants, genes_of_interest){
  
  # MAKE SURE THAT THE INPUT VARIANTS ARE JUST TWO COLUMNS: 1) THE VARIANT COLUMN AND 2) THE GENE COLUMN

  
  
  # Update the colnames if the input background and significant variant lists
  colnames(annotated_significant_variants) <- c('Variant', 'Gene')
  colnames(annotated_background_variants) <- c('Variant', 'Gene')
      
  # Modify significant variant list to have 1 gene per SNV per row
  sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(annotated_significant_variants, splitCols = 'Gene', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
        
  # Modify background Variant list to have 1 gene per SNV per row
  background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(annotated_background_variants, splitCols = 'Gene', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
      
  # Calculate proportions in the significant variant list
      
      # Find index overlap between genes linked to sig variants and a geneset of interest
      sig.modified1.overlap <- which(sig.snvs.modified1$Gene %in% genes_of_interest)
      
      # IF THERE IS NO OVERLAP, ASSIGN 0 TO THE AMOUNT OF OVERLAP
      if (length(sig.modified1.overlap) == 0) {
        
          # Assign NULL to the overlap
          sig.modified1.overlap <- NULL
          
          # number of unique variants near geneset of interest
          sig.in.annotation <- 0
          
          # number of unique variants not near geneset of interest
          sig.notin.annotation <- nrow(annotated_significant_variants) - sig.in.annotation
          
          # fraction of variants near geneset of interest
          sig.fraction.annotation <- sig.in.annotation / (sig.in.annotation + sig.notin.annotation)
        
      } else {
        
          # Extract overlap
          sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
          
          # collect genes for shared variants, if there are any
          sig.modified1.overlap <- aggregate(Gene ~ Variant, sig.modified1.overlap, FUN = toString)
          
          # Order results by alphabetical gene name
          sig.modified1.overlap <- sig.modified1.overlap[order(sig.modified1.overlap$Gene), ]
          
          # number of unique variants near geneset of interest
          sig.in.annotation <- nrow(sig.modified1.overlap)
          
          # number of unique variants not near geneset of interest
          sig.notin.annotation <- nrow(annotated_significant_variants) - sig.in.annotation
          
          # fraction of variants near geneset of interest
          sig.fraction.annotation <- sig.in.annotation / (sig.in.annotation + sig.notin.annotation)
          
      }
  
  # Calculate proportions in the background variant list
      
      # Find index overlap between background snvs and a geneset of interest
      background.modified1.overlap <- which(background.snvs.modified1$Gene %in% genes_of_interest)
      
      # Extract overlap
      background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
      
      # collect genes for shared variants, if there are any
      background.modified1.overlap <- aggregate(Gene ~ Variant, background.modified1.overlap, FUN = toString)
      
      # number of unique variants near geneset of interest
      background.in.annotation <- nrow(background.modified1.overlap)
      
      # number of unique variants not near geneset of interest
      background.notin.annotation <- nrow(annotated_background_variants) - background.in.annotation
      
      # fraction of variants near geneset of interest
      background.fraction.annotation <- background.in.annotation / (background.in.annotation + background.notin.annotation)
      
  # Collect fractions in a df to use for plot
  fraction.annotation <- data.frame(Fraction = c(sig.fraction.annotation, background.fraction.annotation), Group = c('Significant', 'Background'))
      
  # Generate a contingency table as a matrix to run stats
  my.contingency <- as.matrix(data.frame(Not_Annotation = c(sig.notin.annotation, background.notin.annotation), Annotation = c(sig.in.annotation, background.in.annotation)))
  
  # Run Fischer's Exact Test
  my.Fischer.res <- fisher.test(my.contingency)
  
  
  
  # Output 1) the df with fraction of variants with annotation, 2) the Fischer test results, 3) and the df with significant variants that overlap annotation
  return(list(fraction.annotation, my.Fischer.res, sig.modified1.overlap))
  
  
} # END FUNCTION
