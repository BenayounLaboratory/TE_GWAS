Individual_to_CMplot <- function(gwas.results){
  
  # Add an index columns to the results
  gwas.results$index <- 1:nrow(gwas.results)
  
  # Subset the GWAS df and define the CMplot object.
  CMplot_df <- gwas.results[, c('index', 'CHR', 'BP', 'P', 'SNP')]
  
  # Return the final object to plot
  return(CMplot_df)
  
} # END FUNCTION\

Meta_to_CMplot <- function(meta.results){
  
  # Add an index columns to the results
  meta.results$index <- 1:nrow(meta.results)
  
  # Subset the GWAS df and define the CMplot object.
  CMplot_df <- meta.results[, c('index', 'CHR', 'BP', 'P(R)', 'SNP')]
  
  # Return the final object to plot
  return(CMplot_df)
  
} # END FUNCTION

Ensembl_to_symbol <- function(organism, Ensembl_genes){
  
  # THIS FUNCTION MAPS AN ENSEMBL GENE NAME TO ITS SYMBOL
  
  
  
  # Define organism specific parameters
  if (organism == 'hs') {
    organism_dataset <- "hsapiens_gene_ensembl"
    Ensembl_version <- 99
    
  } else if (organism == 'mm') {
    organism_dataset <- "mmusculus_gene_ensembl"
    Ensembl_version <- 99 # DOUBLE CHECK
  }
  

  
  # Retrieve ensembl database
  ensembl = biomaRt::useEnsembl(biomart="ensembl", dataset = organism_dataset, host = 'https://www.ensembl.org', version = Ensembl_version, verbose = TRUE)
  
  # Map Ensembl names to gene symbols
  Mapped_gene_info <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "entrezgene_id"), 
                                     filters = c("ensembl_gene_id"), 
                                     values = Ensembl_genes,
                                     mart = ensembl,
                                     verbose = FALSE,
                                     uniqueRows = TRUE)  
  
  # If genes have multiple mappings, keep only the first. 
  Mapped_gene_info <- Mapped_gene_info[!duplicated(Mapped_gene_info$ensembl_gene_id), ]
  
  # Assign Ensembl names to the rownames
  rownames(Mapped_gene_info) <- Mapped_gene_info$ensembl_gene_id
  
  # Remove ensembl gene id column since they're in rownames
  Mapped_gene_info <- Mapped_gene_info[, -c(1)] 
  
  # Return table with the mapping info
  return(Mapped_gene_info)
  
  
} # END FUNCTION

