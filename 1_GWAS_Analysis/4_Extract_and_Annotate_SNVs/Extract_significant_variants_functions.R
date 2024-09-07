empiricalFDR_to_pval <- function(data.file, null_data.dir, null_file_extension, GWAS_type, arbitrary_pval, empirical_FDR_target){
  
  # THIS FUNCTION TAKES A GWAS RESULTS FILE, AND PERMUTATION GWAS RESULTS FILES, AND DETERMINES THE PVALUE AT WHICH THE EMPIRICAL_FDR_TARGET IS MET
  # null_file_extension will be used to scan a directory for the null files of interest
  # GWAS_type should be 'individual' or 'meta', depending on the GWAS results type
  # abitrary pval will be used to trim GWAS data files so that they can all be loaded into memory
  # empirical FDR target is the empirical FDR we will select a pval threshold for
  # average empirical FDR at a given pvalue is defined as (number_permutation_snps/N_of_permutations)/number_real_snps
  # For meta-analysis, pvalue filtering is on the RANDOM EFFECT pvalue, not the fixed effects pvalue
  
  # Ouput for the function below is 
  # 1) first pvalue where empirical FDR < empirical.FDR.threshold 
  # 2) the exact FDR at that pvalue 
  # 3) the average number of permutation datapoints with pvalue <= the pvalue calculated in result 1 ie false positives and 
  # 4) the number of real data points with pvalue <= the pvalue calculated in result 1. Dividing result 3 by result 4 yields the empirical FDR at that pvalue
  # Note, the function starts at 'arbitrary_pval', calculates the empirical FDR using the real and the permutation data, and consecutively goes to the next smallest pval until the target empirical FDR is reached. If it is not reached, it will spit out NULL.
  
  
  
  # Load the regression file of interest with real/non-permutation results
  GWAS.stats <- fread(file = data.file, header = TRUE, sep = " ", data.table = FALSE)
  
  
  if (GWAS_type == c('individual')) {
    
      # Only keep stats passing arbitrary pval threshold (to avoid loading everything into memory)
      GWAS.stats <- GWAS.stats[which(GWAS.stats$P <= arbitrary_pval), ]
      
      # Define permutation data files extension
      permutation_file_extension <- paste('\\', null_file_extension, sep = '')
  } 

  if (GWAS_type == c('meta')) {
    
      # Only keep stats passing arbitrary pval threshold (to avoid loading everything into memory)
      GWAS.stats <- GWAS.stats[which(GWAS.stats$`P(R)` <= arbitrary_pval), ]
      
      # Define permutation data files extension
      permutation_file_extension <- paste('\\', null_file_extension, '$', sep = '')
  }
  

  # Define permutation data files
  permutations.files <- list.files(null_data.dir, permutation_file_extension, recursive = TRUE, full.names = TRUE)
  
  # Initiate counter 
  permutation_counter <- 0
  
  # Aggregate stats passing the arbitrary pval threshold across ALL permutation files
  for (ith_permutation_file in permutations.files) {
   
    # Update permutation counter
    permutation_counter <- permutation_counter + 1
    
    # Read external file 
    null.GWAS.stats <- fread(file = ith_permutation_file, header = TRUE, sep = " ", data.table = FALSE)
    
    # Only keep stats (passing arbitrary pval threshold)
    if (GWAS_type == c('individual')) {
      
        null.GWAS.stats <- null.GWAS.stats[which(null.GWAS.stats$P <= arbitrary_pval), ]
    } 

    if (GWAS_type == c('meta')) {
      
        null.GWAS.stats <- null.GWAS.stats[which(null.GWAS.stats$`P(R)` <= arbitrary_pval), ]
    }
    
    # If it's the first permutation, define a new dataframe. Otherwise, append new permutation to the growing dataframe.
    if (permutation_counter == 1) {
      
      all_permutation_data <- null.GWAS.stats
      
    } else {
      
      all_permutation_data <- rbind(all_permutation_data, null.GWAS.stats)
    }
    
     
  } # End for loop over permutation files
  

  
  # Order **ALL DATA* by ascending pvalue. NOTE: here, the output from either individual superpopulation or meta-analysis results will use the same variables, to keep downstream code non-redundant
  if (GWAS_type == c('individual')) {
    
      # Order permutation data
      all_permutation_data <- all_permutation_data[order(all_permutation_data$P, decreasing = FALSE), ]
      
      # Order real data
      GWAS.stats <- GWAS.stats[order(GWAS.stats$P, decreasing = FALSE), ]
  } 
  
  if (GWAS_type == c('meta')) {
    
      # Order permutation data
      all_permutation_data <- all_permutation_data[order(all_permutation_data$`P(R)`, decreasing = FALSE), ]
      
      # Order real data
      GWAS.stats <- GWAS.stats[order(GWAS.stats$`P(R)`, decreasing = FALSE), ]
  }
  
  
  
  
  # Start testing pvalues, going from highest to lowest (will keep the first pvalue where FDR < empirical_FDR_target)
  
  # Define the total number of real pvalues to test
  total_pvalues <- nrow(GWAS.stats)
  
  # Calculate the number of permutations
  N_of_permutations <- length(permutations.files)
      
  # Loop over pvalues, from highest to smallest, to home in on the pvalue where FDR is less than the desired target
  for (ith_entry in (total_pvalues:1)) {
    
      if (GWAS_type == c('individual')) {
        
            # Extract the ith pvalue
            ith_pvalue <- GWAS.stats[ith_entry, 'P']
            
            # Count the number of snps with equal or smaller pvalue in the *real* data
            number_real_snps <- nrow(GWAS.stats[which(GWAS.stats$P <= ith_pvalue), ])
            
            # Count the number of snps with equal or smaller pvalue in the *permutation* data
            number_permutation_snps <- nrow(all_permutation_data[which(all_permutation_data$P <= ith_pvalue), ])
      } 
    
      if (GWAS_type == c('meta')) {
        
            # Extract the ith pvalue
            ith_pvalue <- GWAS.stats[ith_entry, 'P(R)']
            
            # Count the number of snps with equal or smaller pvalue in the *real* data
            number_real_snps <- nrow(GWAS.stats[which(GWAS.stats$`P(R)` <= ith_pvalue), ])
            
            # Count the number of snps with equal or smaller pvalue in the *permutation* data
            number_permutation_snps <- nrow(all_permutation_data[which(all_permutation_data$`P(R)` <= ith_pvalue), ])
      }
      
      # Obtain the average number of permutation snps (by dividing by the number of permutations)
      avg_number_permutation_snps <- number_permutation_snps/N_of_permutations
      
      # Calculate empirical FDR at the ith_pvalue
      ith_empirical_FDR <- avg_number_permutation_snps/number_real_snps
      
      # Provide an output if the target FDR is reached
      if (ith_empirical_FDR < empirical_FDR_target) {
        
        # Define the output vector
        output_vector <- c(ith_pvalue, ith_empirical_FDR, avg_number_permutation_snps, number_real_snps)
        
        return(output_vector)
        
      }
    
  } # End for loop over pvalues
  
  
} # END FUNCTION

fill_in_individual_ethnic_group_pval <- function(input_dataframe, ethnic_group_GWAS_file_path, name_of_column_to_fill) {
  
  # FUNCTION INFO: Takes the GWAS results from an individual ethnic group and adds the pvalue to a dataframe containing meta-analysis results
  # input_dataframe is the dataframe to be filled in. 
  
  
  
  
  # Assign RsID to the input_dataframe rownames (for easier manipulation)
  rownames(input_dataframe) <- input_dataframe$SNP
  
  # Load the GWAS results file for an individual superpopulation
  individual.GWAS.stats <- read.csv(file = ethnic_group_GWAS_file_path, header = TRUE, sep = "")
  
      # Keep only individual snps that are in the meta-analysis results (ie are in the input dataframe; it's useful to trim the individual GWAS results early on since the file size is huge)
      individual.GWAS.stats <- individual.GWAS.stats[which(individual.GWAS.stats$SNP %in% input_dataframe$SNP), ]
      
      # Remove duplicate RsIDs from the invidiual GWAS results (if they exist)
      individual.GWAS.stats <- individual.GWAS.stats[!duplicated(individual.GWAS.stats$SNP), ]
      
  # Add pvalues to the input_dataframe
  input_dataframe[individual.GWAS.stats$SNP, name_of_column_to_fill] <- individual.GWAS.stats[, 'P']
  
  # Replace input dataframe rownames with row indices
  rownames(input_dataframe) <- 1:nrow(input_dataframe)


  # Output the final dataframe
  return(input_dataframe)
  
  
} # END FUNCTION

Extract_sig_snps_from_GWAS <- function(output_location, GWAS_stats_df, vector_of_significance_thresholds, output.file.label){
  
  # THIS FUNCTION TAKES ****PVALUES**** CORRESPONDING TO MULTIPLE SIGNFICANT THRESHOLDS (BH FDR, EMPIRICAL FDR, GENOMEWIDE THRESHOLD, ETC.), CHOOSES THE STRICTEST ONE, EXTRACTS SNPS PASSING THAT THRESHOLD FROM A DF WITH META-ANALYSIS RANDOM EFFECT MODEL RESULTS, AND SAVES THEM (IF THERE IS AT LEAST ONE)
  # THIS SCRIPT CURRENTLY CANNOT BE USED FOR INDIVIDUAL SUPER POPULATION GWAS RESULTS, ONLY THE META ANALYSIS
  # THIS FUNCTION WILL ALSO ADD HG19 COORDINATES TO SIGNIFICANT SNVS, SINCE THESE MIGHT BE NECESSARY FOR COLLABORATORS (requires the biomart library)
  
  

  # FILTER VARIANTS
  
  # Define the strictest significance threshold
  final.significance.threshold <- min(vector_of_significance_thresholds)
  
  # Order GWAS results by pvalue (from smallest to largest)
  GWAS_stats_df <- GWAS_stats_df[order(GWAS_stats_df$`P(R)`, decreasing = FALSE), ]
  
  # Extract snps passing FDR (note: <= is used here because the pvalue thresholds are already at FDR < 0.05)
  GWAS.stats.sig <- GWAS_stats_df[GWAS_stats_df$`P(R)` <= final.significance.threshold, ]
  
  
  
  
  
  # MAP TO HG19
  
  # Add columns to hold coorindates
  GWAS.stats.sig$HG19_CHR <- NA
  GWAS.stats.sig$HG19_BP <- NA
  
  # Obtain hg19 coordinates for GWAS RsIDs (NOTE: variants without an RsID will not map)
    
    # Select mart
    ensembl.snp.mart <- useEnsembl("snp", dataset = "hsapiens_snp", GRCh = "37", verbose = TRUE)
    
    # get genomic position
    snp_hg19_mapping_table <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"), 
                                    filters = "snp_filter", 
                                    values = GWAS.stats.sig$SNP, 
                                    mart = ensembl.snp.mart, 
                                    uniqueRows = TRUE,
                                    verbose = FALSE)

    # Define mappings that don't fall on PATCH positions and subset those
    nonpatch.rows <- which(!grepl('PATCH', snp_hg19_mapping_table$chr_name))
    snp_hg19_mapping_table <- snp_hg19_mapping_table[nonpatch.rows, ]

    # order the mapping table by chromosome number (this is done to move patch positions to the end of the list, and preferentially remove them if a non-patched position exists)
    snp_hg19_mapping_table <- snp_hg19_mapping_table[order(snp_hg19_mapping_table$chr_name, decreasing = FALSE), ]
    
    # Find chromosome names of duplicate RsIds (they should correspond to genome patches)
    duplicate_rsid_chromosomes <- unique(snp_hg19_mapping_table[duplicated(snp_hg19_mapping_table$refsnp_id), 'chr_name'])
    
    # remove duplicate RsId entries 
    snp_hg19_mapping_table <- snp_hg19_mapping_table[!duplicated(snp_hg19_mapping_table$refsnp_id), ]
    
# Add hg19 coordinates to GWAS meta-analysis results
    
    # Assign RsIDs as rownames for easier manipulation of the snp_info table
    rownames(GWAS.stats.sig) <- GWAS.stats.sig$SNP
    
    # Add the hg19 chromosome info to the table
    GWAS.stats.sig[snp_hg19_mapping_table$refsnp_id, 'HG19_CHR'] <- snp_hg19_mapping_table$chr_name
    
    # Add the hg19 base position info to the table
    GWAS.stats.sig[snp_hg19_mapping_table$refsnp_id, 'HG19_BP'] <- snp_hg19_mapping_table$chrom_start

  
  
  


    
    
    
    
  # SAVE RESULTS
  
  # Define unique, significant RsIDs
  unique.sig.RsIDs <- unique(GWAS.stats.sig$SNP)
  
  # Define positions for unique, significant RsIDs
  unique.variant.positions <- GWAS.stats.sig[, c('CHR', 'BP')]
  
  # Save, if there are significant snps
  if (nrow(GWAS.stats.sig) > 0) {
    
    # Save statistics for significant SNVs
    write.table(GWAS.stats.sig, file = paste(output_location, output.file.label, '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    # Save unique RsID list
    write.table(unique.sig.RsIDs, file = paste(output_location, output.file.label, '_Unique_RsIDs', '.txt', sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)
    
    # Save unique variant positions
    write.table(unique.variant.positions, file = paste(output_location, output.file.label, '_Unique_Variant_Positions', '.txt', sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)
      
  }
  
  
  
  
  
  # End function
  return(NULL)
  
  
} # END FUNCTION
