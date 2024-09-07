Ensembl_to_symbol <- function(organism, Ensembl_genes){
  
  # needs the following libraries: biomaRt
  
  
  
  # Define organism specific parameters
  if (organism == 'hs') {
    organism_dataset <- "hsapiens_gene_ensembl"
    Ensembl_version <- 99
    
  } else if (organism == 'mm') {
    organism_dataset <- "mmusculus_gene_ensembl"
    Ensembl_version <- 99 # DOUBLE CHECK
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

Extract_DESeq_stats <- function(DESeq_results, padj_limit, organism, output.dir, output_file_prefix_all, output_file_prefix_sig){
  
  # FUNCTION NOTE: Extract statistics from DESeq Results object, remove genes with NA for padj, save padj filtered and unfiltered list. The unfiltered list can be used as gene background.
  
  
  
  # Check if output directory exists. If not, make it. 
  if ( dir.exists(output.dir) == FALSE) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  
  
  # Extract statistics table from DESeq results
  DESeq_stats <- as.data.frame(DESeq_results@listData)
  
  # Assign rownames (genes) to statistics table
  rownames(DESeq_stats) <- DESeq_results@rownames
  
  # Remove rows with NAs
  DESeq_stats <- na.omit(DESeq_stats)
  
  # Getting alternative gene names for remaining genes
  alternative_names <- Ensembl_to_symbol(organism = organism, Ensembl_genes = rownames(DESeq_stats))
  
  # Make columns to hold alternative gene names
  DESeq_stats$external_gene_name <- NA 
  DESeq_stats$hgnc_symbol <- NA
  DESeq_stats$entrezgene_id <- NA
  
  # Assign alternative names
  DESeq_stats[rownames(alternative_names),c('external_gene_name', 'hgnc_symbol', 'entrezgene_id')] <- alternative_names[,]

  # Filter out DEGs
  DESeq_stats.sig <- DESeq_stats[DESeq_stats$padj < padj_limit, ]
  
  # Save stats for the full and significant gene names
  write.table(DESeq_stats.sig, file = paste(output.dir, output_file_prefix_sig, '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)
  write.table(DESeq_stats, file = paste(output.dir, output_file_prefix_all, '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

  
  
  # Output the significant gene stats
  return(DESeq_stats.sig)
  
} # END FUNCTION

Run_DESeq <- function(all_covariates, variant_ID, target_population, counts_df, my.output.dir, my.output.prefix){
  
  
  # This function runs DESeq comparing the genotypes of a variant (SNV OR SV) of interest in either EUR or YRI samples from the Geuvadis project
  # 'population' should be either 'EUR' or 'YRI'
  
  
  
  
  # Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
  padj_limit <- 0.05 
  
  
  if (target_population == 'YRI') {
    
      # Define the final set of covariates to include in the DESeq model
      SampleInfo <- data.frame(row.names = rownames(all_covariates),
                               sex = as.factor(all_covariates$sex),
                               lab = as.factor(all_covariates$lab),
                               PC1 = as.numeric(all_covariates$PC1),
                               PC2 = as.numeric(all_covariates$PC2),
                               EBV = as.numeric(all_covariates$EBV_expr_VST),
                               variant = as.numeric(all_covariates[, variant_ID])
                               )
      
      # Create DESeq2 object
      dds <- DESeqDataSetFromMatrix(countData = counts_df,
                                    colData = SampleInfo,
                                    design = ~ sex + lab + PC1 + PC2 + EBV + variant) 
      
  } else if (target_population == 'EUR') {
    
      # Define the final set of covariates to include in the DESeq model
      SampleInfo <- data.frame(row.names = rownames(all_covariates),
                               sex = as.factor(all_covariates$sex),
                               lab = as.factor(all_covariates$lab),
                               ancestry = as.factor(all_covariates$ancestry),
                               PC1 = as.numeric(all_covariates$PC1),
                               PC2 = as.numeric(all_covariates$PC2),
                               EBV = as.numeric(all_covariates$EBV_expr_VST),
                               variant = as.numeric(all_covariates[, variant_ID])
                               )
      
      # Create DESeq2 object
      dds <- DESeqDataSetFromMatrix(countData = counts_df,
                                    colData = SampleInfo,
                                    design = ~ sex + lab + ancestry + PC1 + PC2 + EBV + variant) 
      
    }
  
  
  # run DESeq2
  dds <- DESeq(dds, parallel = TRUE)

  # Extract DESeq results. 
  DESeq.res <- results(dds, name = c("variant"), alpha = padj_limit, independentFiltering = TRUE)

  # DESeq Stats
  summary(DESeq.res, alpha = padj_limit)

  # Extract significant results and save.
  DESeq.res.sig <- Extract_DESeq_stats(DESeq_results = DESeq.res,
                                       padj_limit = padj_limit,
                                       organism = 'hs',
                                       output.dir = my.output.dir,
                                       output_file_prefix_all = paste(my.output.prefix, '_All_Genes', sep = ''),
                                       output_file_prefix_sig = paste(my.output.prefix, '_FDR5', sep = ''))
  

  
}# END FUNCTION