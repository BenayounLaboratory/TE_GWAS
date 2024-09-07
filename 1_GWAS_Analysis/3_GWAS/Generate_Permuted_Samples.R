# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/Generate_Permuted_Samples_Functions.R")

# Load libraries
library(arrangements)



# PREPARE PERMUTED PHENOTYPES TO GENERATE A NULL DISTRIBUTION FOR GWAS --------------------------------------------


# Define output directories
dir.output.AFR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/AFR/'
dir.output.AMR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/AMR/'
dir.output.EAS <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/EAS/'
dir.output.EUR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/EUR/'
dir.output.SAS <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/SAS/'

# Load sample meta data
sample.phenotypes <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Stratify samples by population
sample.phenotypes.AFR <- sample.phenotypes[which(sample.phenotypes$Super_Population == c('AFR')), ]
sample.phenotypes.AMR <- sample.phenotypes[which(sample.phenotypes$Super_Population == c('AMR')), ]
sample.phenotypes.EAS <- sample.phenotypes[which(sample.phenotypes$Super_Population == c('EAS')), ]
sample.phenotypes.EUR <- sample.phenotypes[which(sample.phenotypes$Super_Population == c('EUR')), ]
sample.phenotypes.SAS <- sample.phenotypes[which(sample.phenotypes$Super_Population == c('SAS')), ]
    
# Permute sample names
permutation.AFR <- generate_permutations(vector_to_permute = rownames(sample.phenotypes.AFR),
                                                              unique_only = T,
                                                              perm_length = nrow(sample.phenotypes.AFR),
                                                              max_permutations = 1e5,
                                                              perm_layout = 'column',
                                                              perms_to_keep = 20) 

permutation.AMR <- generate_permutations(vector_to_permute = rownames(sample.phenotypes.AMR),
                                                              unique_only = T,
                                                              perm_length = nrow(sample.phenotypes.AMR),
                                                              max_permutations = 1e5,
                                                              perm_layout = 'column',
                                                              perms_to_keep = 20) 

permutation.EAS <- generate_permutations(vector_to_permute = rownames(sample.phenotypes.EAS),
                                                              unique_only = T,
                                                              perm_length = nrow(sample.phenotypes.EAS),
                                                              max_permutations = 1e5,
                                                              perm_layout = 'column',
                                                              perms_to_keep = 20) 

permutation.EUR <- generate_permutations(vector_to_permute = rownames(sample.phenotypes.EUR),
                                                              unique_only = T,
                                                              perm_length = nrow(sample.phenotypes.EUR),
                                                              max_permutations = 1e5,
                                                              perm_layout = 'column',
                                                              perms_to_keep = 20)   

permutation.SAS <- generate_permutations(vector_to_permute = rownames(sample.phenotypes.SAS),
                                                              unique_only = T,
                                                              perm_length = nrow(sample.phenotypes.SAS),
                                                              max_permutations = 1e5,
                                                              perm_layout = 'column',
                                                              perms_to_keep = 20) 
    
    
# Generate Plink phenotype permutation files

# AFR
permute_GWAS_phenotypes_for_Plink(phenotypes_df = sample.phenotypes.AFR, 
                                  permutations_df = permutation.AFR, 
                                  phenotype_colname = c('Combined_Plink_Binary'), 
                                  output.dir = dir.output.AFR)

# AMR
permute_GWAS_phenotypes_for_Plink(phenotypes_df = sample.phenotypes.AMR, 
                                  permutations_df = permutation.AMR, 
                                  phenotype_colname = c('Combined_Plink_Binary'), 
                                  output.dir = dir.output.AMR)

 # EAS
permute_GWAS_phenotypes_for_Plink(phenotypes_df = sample.phenotypes.EAS, 
                                  permutations_df = permutation.EAS, 
                                  phenotype_colname = c('Combined_Plink_Binary'), 
                                  output.dir = dir.output.EAS)

# EUR
permute_GWAS_phenotypes_for_Plink(phenotypes_df = sample.phenotypes.EUR, 
                                  permutations_df = permutation.EUR, 
                                  phenotype_colname = c('Combined_Plink_Binary'), 
                                  output.dir = dir.output.EUR)

# SAS
permute_GWAS_phenotypes_for_Plink(phenotypes_df = sample.phenotypes.SAS, 
                                  permutations_df = permutation.SAS, 
                                  phenotype_colname = c('Combined_Plink_Binary'), 
                                  output.dir = dir.output.SAS)
    
        

# Clean the environment
rm(list=ls()) 

