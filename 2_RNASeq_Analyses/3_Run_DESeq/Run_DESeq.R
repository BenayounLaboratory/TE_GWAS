# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/3_Run_DESeq/Run_DESeq_functions.R") 

# Load libraries
library(biomaRt) # to map Ensembl to Entrez
library(DESeq2) # For differential expression 
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(BEDMatrix) # to load genotypes


# Define output directories
my.output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/3_Run_DESeq/DESeq_Results/'  

  

# CASES VS CONTROLS ------------------------------------------------------------------------------------------------------------------------------
  

# LOAD COUNTS AND KNOWN COVARIATES

# Load filtered counts (USE COUNTS FROM MY EQTL STUDY)
counts.YRI <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/All_counts_YRI_86_filtered.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
counts.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/All_counts_EUR_358_filtered.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load the covariates data
covariates_table.YRI  <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/COVARIATES_YRI.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.YRI <- covariates_table.YRI[colnames(counts.YRI), ]      

covariates_table.EUR  <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/COVARIATES_EUR.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
covariates_table.EUR <- covariates_table.EUR[colnames(counts.EUR), ]   

# Load case/control status
singleton.counts <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Only keep samples with transcriptomes
    singleton.counts.EUR <- singleton.counts[which(rownames(singleton.counts) %in% colnames(counts.EUR)), ]
    singleton.counts.YRI <- singleton.counts[which(rownames(singleton.counts) %in% colnames(counts.YRI)), ]
    
    # Organize singleton counts in the same order as counts and covariates
    singleton.counts.YRI <- singleton.counts.YRI[rownames(covariates_table.YRI), ]
    singleton.counts.EUR <- singleton.counts.EUR[rownames(covariates_table.EUR), ]
    
    # Add a column to hold case/control group
    singleton.counts.EUR$GWAS_Group <- 'control'
    singleton.counts.YRI$GWAS_Group <- 'control'
    
    # Update the cases
    singleton.counts.EUR[which(singleton.counts.EUR$Combined_Plink_Binary == 2), 'GWAS_Group'] <- 'case'
    singleton.counts.YRI[which(singleton.counts.YRI$Combined_Plink_Binary == 2), 'GWAS_Group'] <- 'case'

# Add case/control covariates table
covariates_table.EUR$Combined_singletons <- singleton.counts.EUR$GWAS_Group
covariates_table.YRI$Combined_singletons <- singleton.counts.YRI$GWAS_Group

# Add global singleton number to the covariates table
covariates_table.EUR$Combined_singletons_count <- as.numeric(singleton.counts.EUR$Combined_unique_insertions)
covariates_table.YRI$Combined_singletons_count <- as.numeric(singleton.counts.YRI$Combined_unique_insertions)  
  
    



# RUN DESEQ

# Define alpha: the significance cutoff used for optimizing the independent filtering (by default 0.1).   If the adjusted p-value cutoff (FDR) will be a value other than 0.1, alpha should be set to that value.
padj_limit <- 0.05 

# Define the final set of covariates to include in the DESeq model
SampleInfo.EUR <- data.frame(row.names = rownames(covariates_table.EUR),
                             sex = as.factor(covariates_table.EUR$sex),
                             lab = as.factor(covariates_table.EUR$lab),
                             ancestry = as.factor(covariates_table.EUR$ancestry),
                             PC1 = as.numeric(covariates_table.EUR$PC1),
                             PC2 = as.numeric(covariates_table.EUR$PC2),
                             EBV = as.numeric(covariates_table.EUR$EBV_expr_VST),
                             Combined_singletons = as.factor(covariates_table.EUR$Combined_singletons)
                             )

SampleInfo.YRI <- data.frame(row.names = rownames(covariates_table.YRI),
                             sex = as.factor(covariates_table.YRI$sex),
                             lab = as.factor(covariates_table.YRI$lab),
                             PC1 = as.numeric(covariates_table.YRI$PC1),
                             PC2 = as.numeric(covariates_table.YRI$PC2),
                             EBV = as.numeric(covariates_table.YRI$EBV_expr_VST),
                             Combined_singletons = as.factor(covariates_table.YRI$Combined_singletons)
                             )

# Create DESeq2 object
dds.EUR <- DESeqDataSetFromMatrix(countData = counts.EUR,
                                  colData = SampleInfo.EUR,
                                  design = ~ sex + lab + ancestry + PC1 + PC2 + EBV + Combined_singletons) 

dds.YRI <- DESeqDataSetFromMatrix(countData = counts.YRI,
                                  colData = SampleInfo.YRI,
                                  design = ~ sex + lab + PC1 + PC2 + EBV + Combined_singletons) 

# run DESeq2
dds.EUR <- DESeq(dds.EUR, parallel = TRUE)
dds.YRI <- DESeq(dds.YRI, parallel = TRUE)

# Extract DESeq results. 
res.EUR <- results(dds.EUR, contrast = c("Combined_singletons", "case", "control"), alpha = padj_limit, independentFiltering = TRUE)
res.YRI <- results(dds.YRI, contrast = c("Combined_singletons", "case", "control"), alpha = padj_limit, independentFiltering = TRUE)

# DESeq Stats
summary(res.EUR, alpha = padj_limit)
summary(res.YRI, alpha = padj_limit)

# Extract significant results and save.
res.EUR.sig <- Extract_DESeq_stats(DESeq_results = res.EUR, 
                                   padj_limit = padj_limit, 
                                   organism = 'hs', 
                                   output.dir = my.output.dir,
                                   output_file_prefix_all = paste('Cases_vs_Controls_EUR_All_Genes', sep = ''),
                                   output_file_prefix_sig = paste('Cases_vs_Controls_EUR_FDR5', sep = ''))

res.YRI.sig <- Extract_DESeq_stats(DESeq_results = res.YRI, 
                                   padj_limit = padj_limit, 
                                   organism = 'hs', 
                                   output.dir = my.output.dir,
                                   output_file_prefix_all = paste('Cases_vs_Controls_YRI_All_Genes', sep = ''),
                                   output_file_prefix_sig = paste('Cases_vs_Controls_YRI_FDR5', sep = ''))



# POLYMORPHIC REPEATS ------------------------------------------------------------------------------------------------------------------------------


# DEFINE GENOTYPES FOR TOP POLYMORPHIC REPEATS

# Specify the genotype BED files
plink_file_path.YRI <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/Final_SNVs_SVs_AFR.bed'
plink_file_path.EUR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes/Final_SNVs_SVs_EUR.bed'

# Make R object to stream snp genotypes into memory (rows are samples and columns are snps). # NOTE: If used, MT snps will need to be assigned a snpid, since they dont have rsid.
binary_genotypes.YRI <- BEDMatrix(path = plink_file_path.YRI, simple_names = T)
binary_genotypes.EUR <- BEDMatrix(path = plink_file_path.EUR, simple_names = T)

# Define the snps of interest
top.snps.YRI <- c('INV_delly_INV00066128', 'INV_delly_INV00003623', 'ALU_umary_ALU_7919', 'ALU_umary_ALU_3176', 'ALU_umary_ALU_8971', 'L1_umary_LINE1_1066') 
top.snps.EUR <- c('INV_delly_INV00066128', 'INV_delly_INV00003623', 'ALU_umary_ALU_7919', 'ALU_umary_ALU_3176', 'ALU_umary_ALU_8971', 'L1_umary_LINE1_1066') 

# Extract genotypes for SVs to be tested
top.genotypes.YRI <- as.data.frame(binary_genotypes.YRI[rownames(covariates_table.YRI), top.snps.YRI])
top.genotypes.EUR <- as.data.frame(binary_genotypes.EUR[rownames(covariates_table.EUR), top.snps.EUR])

# Add SVs to covariates table
covariates_table.YRI <- cbind(covariates_table.YRI, top.genotypes.YRI)
covariates_table.EUR <- cbind(covariates_table.EUR, top.genotypes.EUR)


    


# RUN DESEQ

# YRI SVs
Run_DESeq(all_covariates = covariates_table.YRI, 
          variant_ID = c('INV_delly_INV00066128'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('YRI_INV_Number01_INV_delly_INV00066128'))  

Run_DESeq(all_covariates = covariates_table.YRI, 
          variant_ID = c('INV_delly_INV00003623'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('YRI_INV_Number02_INV_delly_INV00003623'))  

Run_DESeq(all_covariates = covariates_table.YRI, 
          variant_ID = c('ALU_umary_ALU_7919'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('YRI_Alu_Number01_ALU_umary_ALU_7919'))  

Run_DESeq(all_covariates = covariates_table.YRI, 
          variant_ID = c('ALU_umary_ALU_3176'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('YRI_Alu_Number02_ALU_umary_ALU_3176'))  

Run_DESeq(all_covariates = covariates_table.YRI, 
          variant_ID = c('ALU_umary_ALU_8971'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('YRI_Alu_Number03_ALU_umary_ALU_8971'))  

Run_DESeq(all_covariates = covariates_table.YRI, 
          variant_ID = c('L1_umary_LINE1_1066'), 
          target_population = c('YRI'), 
          counts_df = counts.YRI, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('YRI_L1_Number01_L1_umary_LINE1_1066'))  


# EUR SVs
Run_DESeq(all_covariates = covariates_table.EUR, 
          variant_ID = c('INV_delly_INV00066128'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('EUR_INV_Number01_INV_delly_INV00066128'))  

Run_DESeq(all_covariates = covariates_table.EUR, 
          variant_ID = c('INV_delly_INV00003623'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('EUR_INV_Number02_INV_delly_INV00003623'))  

Run_DESeq(all_covariates = covariates_table.EUR, 
          variant_ID = c('ALU_umary_ALU_7919'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('EUR_Alu_Number01_ALU_umary_ALU_7919'))  

Run_DESeq(all_covariates = covariates_table.EUR, 
          variant_ID = c('ALU_umary_ALU_3176'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('EUR_Alu_Number02_ALU_umary_ALU_3176'))  

Run_DESeq(all_covariates = covariates_table.EUR, 
          variant_ID = c('ALU_umary_ALU_8971'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('EUR_Alu_Number03_ALU_umary_ALU_8971'))  

Run_DESeq(all_covariates = covariates_table.EUR, 
          variant_ID = c('L1_umary_LINE1_1066'), 
          target_population = c('EUR'), 
          counts_df = counts.EUR, 
          my.output.dir = my.output.dir, 
          my.output.prefix = c('EUR_L1_Number01_L1_umary_LINE1_1066'))  



# SESSION INFO ------------------------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/3_Run_DESeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Run_DESeq.txt", sep =""))
sessionInfo()
sink()      
    

# Clean the environment
rm(list=ls())

