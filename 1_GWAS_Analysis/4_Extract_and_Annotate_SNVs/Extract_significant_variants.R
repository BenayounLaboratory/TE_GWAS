# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(data.table) # for fread and fsave, used in the functions
library(biomaRt) # for mapping RsIds to hg19 positions
library(BEDMatrix) # To stream genotypes for histograms

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Extract_significant_variants_functions.R") 

# Define output directory for significant snp files
output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/'



# DETERMINE EMPIRICAL FDR/PVALUE THRESHOLDS ----------------------------------------------------------------------------------------


# DEFINE INPUT PARAMETERS

# Define the target empirical FDR threshold
empirical.FDR.threshold <- 0.05

# Make a dataframe to hold the pvalue thresholds (and other stats) for the desired empirical FDR threshold
all.empirical.data <- as.data.frame(matrix(nrow = 6, ncol = 4, NA))

# Assign colnames
colnames(all.empirical.data) <- c('P_threshold', 'Empirical_FDR_at_P', 'Avg_Num_Permutation_Points_Below_P', 'Num_Real_Points_Below_P')

# Assign rownames
rownames(all.empirical.data) <- c('Combined_meta', 'Combined_AFR', 'Combined_AMR', 'Combined_EAS', 'Combined_EUR', 'Combined_SAS')





# CALCULATE EMPIRICAL FDR THRESHOLDS 
  
# Combined Singletons, Meta-analysis
empirical_summary_stats <-  empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/plink.meta', 
                                                 null_data.dir = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/',  
                                                 null_file_extension = '.meta',  
                                                 GWAS_type = c('meta'), 
                                                 arbitrary_pval = 1.8e-05, 
                                                 empirical_FDR_target = empirical.FDR.threshold)

    # Print out the results
    empirical_summary_stats # 1.3960e-05 4.9106e-02 3.8450e+01 7.8300e+02
    
    # Add results to the df
    if (!is.null(empirical_summary_stats)) {
        all.empirical.data['Combined_meta', ] <- empirical_summary_stats
    }
      

# Combined Singletons, AFR only
empirical_summary_stats <-  empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AFR_GWAS.assoc.logistic', 
                                                 null_data.dir = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations',  
                                                 null_file_extension = paste('Singletons_AFR_', '\\d+' ,'.assoc.logistic', sep = ''),  
                                                 GWAS_type = c('individual'), 
                                                 arbitrary_pval = 9e-6, 
                                                 empirical_FDR_target = empirical.FDR.threshold)

    # Print out the results
    empirical_summary_stats # 4.606000e-06 4.888889e-02 8.800000e+00 1.800000e+02
    
    # Add results to the df
    if (!is.null(empirical_summary_stats)) {
        all.empirical.data['Combined_AFR', ] <- empirical_summary_stats
    }

    
# Combined Singletons, AMR only
empirical_summary_stats <-  empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AMR_GWAS.assoc.logistic', 
                                                 null_data.dir = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations',  
                                                 null_file_extension = paste('Singletons_AMR_', '\\d+' ,'.assoc.logistic', sep = ''),  
                                                 GWAS_type = c('individual'), 
                                                 arbitrary_pval = 5e-5, 
                                                 empirical_FDR_target = empirical.FDR.threshold)

    # Print out the results
    empirical_summary_stats # 1.069e-06 4.400e-02 2.200e+00 5.000e+01
    
    # Add results to the df
    if (!is.null(empirical_summary_stats)) {
        all.empirical.data['Combined_AMR', ] <- empirical_summary_stats
    }
    
    
# Combined Singletons, EAS only
empirical_summary_stats <-  empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EAS_GWAS.assoc.logistic', 
                                                 null_data.dir = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations',  
                                                 null_file_extension = paste('Singletons_EAS_', '\\d+' ,'.assoc.logistic', sep = ''),  
                                                 GWAS_type = c('individual'), 
                                                 arbitrary_pval = 9e-6, 
                                                 empirical_FDR_target = empirical.FDR.threshold)

    # Print out the results
    empirical_summary_stats # NULL
    
    # Add results to the df
    if (!is.null(empirical_summary_stats)) {
        all.empirical.data['Combined_EAS', ] <- empirical_summary_stats
    }
    
    
# Combined Singletons, EUR only
empirical_summary_stats <-  empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EUR_GWAS.assoc.logistic', 
                                                 null_data.dir = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations',  
                                                 null_file_extension = paste('Singletons_EUR_', '\\d+' ,'.assoc.logistic', sep = ''),  
                                                 GWAS_type = c('individual'), 
                                                 arbitrary_pval = 9.0e-6, 
                                                 empirical_FDR_target = empirical.FDR.threshold)

    # Print out the results
    empirical_summary_stats # 9.164e-08 0.000e+00 0.000e+00 1.000e+00
    
    # Add results to the df
    if (!is.null(empirical_summary_stats)) {
        all.empirical.data['Combined_EUR', ] <- empirical_summary_stats
    }

    
# Combined Singletons, SAS only
empirical_summary_stats <-  empiricalFDR_to_pval(data.file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/SAS_GWAS.assoc.logistic', 
                                                 null_data.dir = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations',  
                                                 null_file_extension = paste('Singletons_SAS_', '\\d+' ,'.assoc.logistic', sep = ''),  
                                                 GWAS_type = c('individual'), 
                                                 arbitrary_pval = 9e-6, 
                                                 empirical_FDR_target = empirical.FDR.threshold)

    # Print out the results
    empirical_summary_stats # 3.097000e-08 2.352941e-02 4.000000e-01 1.700000e+01
    
    # Add results to the df
    if (!is.null(empirical_summary_stats)) {
        all.empirical.data['Combined_SAS', ] <- empirical_summary_stats
    }



    
    
# Save the results
write.table(all.empirical.data, file = paste(output.dir, "EMPIRICAL_PVALUE_AND_FDR_THRESHOLDS", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)

    

# EXTRACT SIGNIFICANT SNVs ----------------------------------------------------------------------------------------


# DEFINE INPUT PARAMETERS

# Load table with empirical pval/FDR thresholds
empirical.thresholds <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/EMPIRICAL_PVALUE_AND_FDR_THRESHOLDS.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Define the desired BH FDR threshold (I will use both empirical and BH FDRs, considering significant those SNVs with pvalues that are smaller than the pval at the stricter threshold)
BH.FDR.threshold <- 0.05

  


  





# ADD PVALUES FROM INDIVIDUAL POPULATIONS

# Read meta-analysis results
meta.results <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/plink.meta', header = TRUE, sep = " ", data.table = FALSE)
      
    # Add BH FDR columns
    meta.results$FDR <- p.adjust(meta.results$`P(R)`, method = 'BH')

# Load table with details of each SNP used in the GWAS analysis
my.SNP.metadata <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_Variant_Map_AFR_SNVs_and_SVs.txt', header = TRUE, sep = "")

    # Only keep metadata for variants in the GWAS meta-analysis
    my.SNP.metadata <- my.SNP.metadata[which(my.SNP.metadata$RsID %in% meta.results$SNP), ]
    
    # Assign RsID as the rowname
    rownames(my.SNP.metadata) <- my.SNP.metadata$RsID

# Since A2 allele info is missing from the GWAS results, fill it in using the REF allele info from my snp metadata table 
# Note, A2 corresponds to the REF allele since I used -keep-allele-order with PLINK
meta.results[, 'A2'] <- my.SNP.metadata[meta.results$SNP, 'REF']

# Add NA columns to fill in the pvalue from each individual superpopulation analysis
meta.results$Pval_AFR = NA
meta.results$Pval_AMR = NA
meta.results$Pval_EAS = NA
meta.results$Pval_EUR = NA
meta.results$Pval_SAS = NA

# Run function to recursively fill-in individual ethnic group pvalues for variants identified in the meta-analysis
    
# AFR
meta.results <- fill_in_individual_ethnic_group_pval(input_dataframe = meta.results,
                                                     ethnic_group_GWAS_file_path = c('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AFR_GWAS.assoc.logistic'),
                                                     name_of_column_to_fill = c('Pval_AFR'))    
# AMR
meta.results <- fill_in_individual_ethnic_group_pval(input_dataframe = meta.results,
                                                     ethnic_group_GWAS_file_path = c('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AMR_GWAS.assoc.logistic'),
                                                     name_of_column_to_fill = c('Pval_AMR'))   
# EAS
meta.results <- fill_in_individual_ethnic_group_pval(input_dataframe = meta.results,
                                                     ethnic_group_GWAS_file_path = c('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EAS_GWAS.assoc.logistic'),
                                                     name_of_column_to_fill = c('Pval_EAS'))   
# EUR
meta.results <- fill_in_individual_ethnic_group_pval(input_dataframe = meta.results,
                                                     ethnic_group_GWAS_file_path = c('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EUR_GWAS.assoc.logistic'),
                                                     name_of_column_to_fill = c('Pval_EUR'))   
# SAS
meta.results <- fill_in_individual_ethnic_group_pval(input_dataframe = meta.results,
                                                     ethnic_group_GWAS_file_path = c('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/SAS_GWAS.assoc.logistic'),
                                                     name_of_column_to_fill = c('Pval_SAS'))   

  


  





# EXTRACT SIGNIFICANT SNPS : Combined Singletons, Meta-analysis

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(meta.results[meta.results$FDR < BH.FDR.threshold, 'P(R)']) 
BH.pval_threshold # 5.995e-06

# Define empirical FDR pvalue
empirical.pval_threshold <- empirical.thresholds[c('Combined_meta'), c('P_threshold')] # 1.396e-05
  
# Run function to extract significant snps and annotate with hg19 coordinates
Extract_sig_snps_from_GWAS(output_location = output.dir, 
                           GWAS_stats_df = meta.results, 
                           vector_of_significance_thresholds = c(BH.pval_threshold, empirical.pval_threshold), 
                           output.file.label = c('Singleton_Meta-Analysis_Significant_Variants'))
            
    

# PLOT ALLELE FREQUENCY HISTOGRAMS FOR SIGNIFICANT SNVs/SVs ----------------------------------------------------------

  
# Load significant GWAS SNVs/SVs
GWAS.sig.SNPs <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/Singleton_Meta-Analysis_Significant_Variants.txt', header = TRUE, sep = "")

# Define vector of ethnic groups to loop over
ethnic_groups <- c('AFR', 'AMR', 'EUR', 'EAS', 'SAS')

# Initiate pdf to hold histograms
pdf(paste(output.dir, "Histogram_Singleton_Significant_Variants_MAFs_per_population", ".pdf", sep=""), width = 4.5, height = 9)

    # Make plot object
    par(mfrow = c(5,1), pty ="m") # 5 rows, 1 columns
    
    # LOOP OVER ETHNIC GROUPS FOR THE SPECIFIED SNP
    for (group_i in ethnic_groups) {
        
          # Define the genotypes to use, changing with each ethnic group
          # NOTE: Rows are human sample IDs, and columns are snps
      
          if (group_i == 'AFR') {
                # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). 
                AFR.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/Final_SNVs_SVs_AFR.bed'
                group_i_genotypes <- BEDMatrix(path = AFR.plink_file_path, simple_names = T)
          }
      
          if (group_i == 'AMR') {
                # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). 
                AMR.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_Genotypes/Final_SNVs_SVs_AMR.bed'
                group_i_genotypes <- BEDMatrix(path = AMR.plink_file_path, simple_names = T)
          }
      
          if (group_i == 'EUR') {
                # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). 
                EUR.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes/Final_SNVs_SVs_EUR.bed'
                group_i_genotypes <- BEDMatrix(path = EUR.plink_file_path, simple_names = T)
          }
      
          if (group_i == 'EAS') {
                # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). 
                EAS.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_Genotypes/Final_SNVs_SVs_EAS.bed'
                group_i_genotypes <- BEDMatrix(path = EAS.plink_file_path, simple_names = T)
          }
      
          if (group_i == 'SAS') {
                # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). 
                SAS.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_Genotypes/Final_SNVs_SVs_SAS.bed'
                group_i_genotypes <- BEDMatrix(path = SAS.plink_file_path, simple_names = T)
          }
      
      
          
      
          # Keep significant snps that are present in the ith ethnic group
          sig_snps_in_ith_population <- intersect(unique(GWAS.sig.SNPs$SNP), colnames(group_i_genotypes))
    
          # Extract significant snp genotypes for the ith population
          ith_pop_sig_snp_genotypes <- group_i_genotypes[, sig_snps_in_ith_population]
          
          # Define the total number of alleles in the ith population
          ith_pop_total_alleles <- 2 * nrow(ith_pop_sig_snp_genotypes)
          
          # Calculate the frequency of the alternate allele (which should be *mostly* the minor allele); defined as (alternate allele count / total alleles)
          ith_pop_sig_snp_MAFs <- colSums(ith_pop_sig_snp_genotypes) / ith_pop_total_alleles
          
                # For snps where the alternate allele is not the minor allele (ie MAF > 0.50), do 1-alternate_frequency to get the actual MAF
          
                    # Define indices where MAF > 0.50
                    snps_where_MAF_above_0.50 <- which(ith_pop_sig_snp_MAFs > 0.50)
                    
                    # Update the MAF at those indices
                    ith_pop_sig_snp_MAFs[snps_where_MAF_above_0.50] <- 1 - ith_pop_sig_snp_MAFs[snps_where_MAF_above_0.50]
          
          # Generate histogram
    
          hist(ith_pop_sig_snp_MAFs,
                main = paste("Significant SNV/SV MAFs in the ", group_i, " Population", sep = ''),
                xlim = c(0, 0.50),
                breaks = seq(0, 0.5, length.out = 101),
                xlab = "MAF",
                ylab = "Number of SNPs",
                col = "blue",
                freq = TRUE # Uses counts
                )
    
          
    } # CLOSE FOR LOOP OVER ETHNIC GROUPS

# Close PDF
dev.off()
            
          
          
# SESSION INFO ----------------------------------------------------------
  

# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Session_Info/'
    
sink(file = paste(dir.session_info, "Session_Info_Extract_Significant_SNVs.txt", sep =""))
sessionInfo()
sink()    

# Clear variables from the environment
rm(list=ls())


