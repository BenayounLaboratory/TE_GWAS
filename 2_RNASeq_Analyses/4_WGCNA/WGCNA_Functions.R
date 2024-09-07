Ensembl_to_symbol <- function(organism, Ensembl_genes){
  
  # needs the following libraries: biomaRt
  
  
  
  # Define organism specific parameters
  if (organism == 'hs') {
    organism_dataset <- "hsapiens_gene_ensembl"
    #Ensembl_version <- 99
    
  } else if (organism == 'mm') {
    organism_dataset <- "mmusculus_gene_ensembl"
    #Ensembl_version <- 99 
  }
  
  
  # Retrieve ensembl database
  #ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset = organism_dataset, host = 'https://www.ensembl.org', version = Ensembl_version, verbose = TRUE)
  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset = organism_dataset, host = 'https://www.ensembl.org', verbose = TRUE)
  
  # Map Ensembl names to gene symbols
  Mapped_gene_info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "entrezgene_id"), 
                                     filters = c("ensembl_gene_id"), 
                                     values = Ensembl_genes,
                                     mart = ensembl,
                                     verbose = FALSE,
                                     uniqueRows = TRUE)  
  
  # If genes have multiple mappings, keep only the first. NOTE: THERE MIGHT BE ISSUES WITH ENTREZGENE_ID MAPPING.
  Mapped_gene_info <- Mapped_gene_info[!duplicated(Mapped_gene_info$ensembl_gene_id), ]
  
  # Cleanup mapping table.
  rownames(Mapped_gene_info) <- Mapped_gene_info$ensembl_gene_id
  Mapped_gene_info <- Mapped_gene_info[, -c(1)] # Remove ensembl gene id column since they're in rownames
  
  # Return table with the mapping info
  return(Mapped_gene_info)
  
  
} # END FUNCTION

Add_GeneSymbols_to_DF <- function(input_df, organism){
  
  # FUNCTION NOTE: Rownames of the input DF should be the EnsemblID
  
  # Getting alternative gene names for genes
  alternative_names <- Ensembl_to_symbol(organism = organism, Ensembl_genes = rownames(input_df))
  
  # Make columns to hold alternative gene names
  input_df$external_gene_name <- NA 
  input_df$hgnc_symbol <- NA
  input_df$entrezgene_id <- NA
  
  # Assign alternative names
  input_df[rownames(alternative_names), c('external_gene_name', 'hgnc_symbol', 'entrezgene_id')] <- alternative_names[, ]
  
  # Output the updated input df
  return(input_df)
  
} # END FUNCTION

trait_cor_meta <- function(consistent.logical.matrix, input.consensus, input.multiset, multiset_type){
  
  # use fisher's method to combine pvalues from the trait-module correlation analysis across consensus networks. also average the correlations for illustration purposes.
  # consistent.logical.matrix is a logical matrix specifying whether a trait-module relationship has the same direction across consensus networks
  # input.consensus is the consensus correlation or pvalue matrix that needs to be filled in
  # input.multiset is the multiset of correlations or pvalues that will be meta-analyzed
  # THIS FUNCTION ASSUMED THAT THERE ARE ONLY 2 ENTRIES IN THE MULTISET (i.e. EUR and YRI)
  
  if (multiset_type == 'correlation') {
      
      # take the mean of correlation values
      input.consensus[consistent.logical.matrix] = (input.multiset[[1]][consistent.logical.matrix] + input.multiset[[2]][consistent.logical.matrix])/2
      
  }
  
  if (multiset_type == 'pvalue') {
    
      # Split multiset into individual sets
      pval.set1 <- input.multiset[[1]][consistent.logical.matrix]
      pval.set2 <- input.multiset[[2]][consistent.logical.matrix]
      
      # Make a vector to hold meta pvalues
      pval.meta <- pval.set1
      
      # Loop over all entries in each individual set
      for (ith_pvalue in 1:length(pval.set1)) {
        
          # Calculate the meta pvalue using Fishers method
          my.meta.pvalue <- poolr::fisher(p = c(pval.set1[ith_pvalue], pval.set2[ith_pvalue]))
          
          # Update the meta pvalue vector
          pval.meta[ith_pvalue] <- my.meta.pvalue$p
        
      }
     
      # Fill in the consensus pvalue matrix
      input.consensus[consistent.logical.matrix] <- pval.meta
    
  }
  

  # output the updated correlations and pvalues
  return(input.consensus)
  
  
} # END FUNCTION
