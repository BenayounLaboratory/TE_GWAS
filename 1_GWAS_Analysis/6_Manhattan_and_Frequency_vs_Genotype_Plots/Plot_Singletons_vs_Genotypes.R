# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(data.table) # for fread and fsave
library(BEDMatrix) # To stream genotypes
library(beeswarm) # for beeswarm plots

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/6_Manhattan_and_Frequency_vs_Genotype_Plots/Plot_Singletons_vs_Genotypes_Functions.R")

# Specify directory for annotated SNVs
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/6_Manhattan_and_Frequency_vs_Genotype_Plots/Singleton_Frequency_vs_Genotype_Plots/'





# LOAD DATA
    
# Load singleton counts per sample
sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load significant annotated GWAS variants (SNVs and SVs)
sig.variants.all <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/All_Annotated_GWAS_Significant_SNVs_AND_SVs.csv', header = TRUE, sep = ",")

    # make a second df to hold non-blacklisted variants
    sig.variants.greenlist <- sig.variants.all
    
    # Remove blacklisted SNPs
    sig.variants.greenlist <- sig.variants.greenlist[which(is.na(sig.variants.greenlist$Blacklist_Label)), ]
    
    # Assign variant ID to rownames
    rownames(sig.variants.greenlist) <- sig.variants.greenlist$SNP

# Load snp mapping info
variant.mapping.info <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_Variant_Map_AFR_SNVs_and_SVs.txt', header = TRUE, sep = "\t", data.table = FALSE)

    # Only keep snps that are in the significant snp table
    variant.mapping.info <- variant.mapping.info[which(variant.mapping.info$RsID %in% sig.variants.greenlist$SNP), ]
    
    # Update rownames
    rownames(variant.mapping.info) <- variant.mapping.info$RsID

# Make R object to stream variant genotypes into memory (rows are samples and columns are snps). 
AFR.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/Final_SNVs_SVs_AFR.bed'
AFR_binary_genotypes <- BEDMatrix(path = AFR.plink_file_path, simple_names = T)

AMR.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_Genotypes/Final_SNVs_SVs_AMR.bed'
AMR_binary_genotypes <- BEDMatrix(path = AMR.plink_file_path, simple_names = T)

EUR.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes/Final_SNVs_SVs_EUR.bed'
EUR_binary_genotypes <- BEDMatrix(path = EUR.plink_file_path, simple_names = T)
   
EAS.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_Genotypes/Final_SNVs_SVs_EAS.bed'
EAS_binary_genotypes <- BEDMatrix(path = EAS.plink_file_path, simple_names = T)

SAS.plink_file_path <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_Genotypes/Final_SNVs_SVs_SAS.bed'
SAS_binary_genotypes <- BEDMatrix(path = SAS.plink_file_path, simple_names = T)





# EXTRACT GREENLIST SIGNIFICANT VARIANT GENOTYPES

# Make dataframe to hold genotypes (rows are samples and columns are snps)
sig.variant.genotypes <- data.frame(matrix(NA, nrow = nrow(sample_info), ncol = nrow(sig.variants.greenlist)))

# Assign colnames
colnames(sig.variant.genotypes) <- sig.variants.greenlist$SNP

# Assign rownames
rownames(sig.variant.genotypes) <- rownames(sample_info)

# Fill in table

# AFR

    # Define indices for significant variants in the bed file
    sig.variant.indices <- which(colnames(sig.variant.genotypes) %in% colnames(AFR_binary_genotypes))
    
    # DefinesIDs for significant variants in the bed file
    sig.variant.ID <- colnames(sig.variant.genotypes[, sig.variant.indices])

    # Extract genotypes for variants in both tables
    sig.variant.genotypes[rownames(AFR_binary_genotypes), sig.variant.ID] <- AFR_binary_genotypes[, sig.variant.ID]
  
# AMR

    # Define indices for significant variants in the bed file
    sig.variant.indices <- which(colnames(sig.variant.genotypes) %in% colnames(AMR_binary_genotypes))
    
    # DefinesIDs for significant variants in the bed file
    sig.variant.ID <- colnames(sig.variant.genotypes[, sig.variant.indices])

    # Extract genotypes for variants in both tables
    sig.variant.genotypes[rownames(AMR_binary_genotypes), sig.variant.ID] <- AMR_binary_genotypes[, sig.variant.ID]
  
# EUR

    # Define indices for significant variants in the bed file
    sig.variant.indices <- which(colnames(sig.variant.genotypes) %in% colnames(EUR_binary_genotypes))
    
    # DefinesIDs for significant variants in the bed file
    sig.variant.ID <- colnames(sig.variant.genotypes[, sig.variant.indices])

    # Extract genotypes for variants in both tables
    sig.variant.genotypes[rownames(EUR_binary_genotypes), sig.variant.ID] <- EUR_binary_genotypes[, sig.variant.ID]    

# EAS

    # Define indices for significant variants in the bed file
    sig.variant.indices <- which(colnames(sig.variant.genotypes) %in% colnames(EAS_binary_genotypes))
    
    # DefinesIDs for significant variants in the bed file
    sig.variant.ID <- colnames(sig.variant.genotypes[, sig.variant.indices])

    # Extract genotypes for variants in both tables
    sig.variant.genotypes[rownames(EAS_binary_genotypes), sig.variant.ID] <- EAS_binary_genotypes[, sig.variant.ID]      

 # SAS

    # Define indices for significant variants in the bed file
    sig.variant.indices <- which(colnames(sig.variant.genotypes) %in% colnames(SAS_binary_genotypes))
    
    # DefinesIDs for significant variants in the bed file
    sig.variant.ID <- colnames(sig.variant.genotypes[, sig.variant.indices])

    # Extract genotypes for variants in both tables
    sig.variant.genotypes[rownames(SAS_binary_genotypes), sig.variant.ID] <- SAS_binary_genotypes[, sig.variant.ID]         

        
# Save file
write.table(sig.variant.genotypes, file = paste(dir.output, "GWAS_Significant_Greenlist_SNV_SV_Genotypes.txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = F)





    
    
    


# PLOT GLOBAL SINGLETON FREQUENCY VS GENOTYPE FOR SIGNIFICANT POLYMORPHIC STRUCTURAL VARIANTS (SVs)
    
# Load significant annotated GWAS variants (SVs only)
sig.SV.greenlist <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/Annotated_Significant_Variants_Greenlist_SVs.csv', header = TRUE, sep = ",")

# Define the variants of interest
interesting_variants <- sig.SV.greenlist$SNP
    
# Extract the info for these snps
interesting_variants.stats <- sig.variants.greenlist[interesting_variants, c('SNP', 'OR.R.', 'FDR')] # manually confirm that genes of interest are linked
    
    # sort results by FDR
    interesting_variants.stats <- interesting_variants.stats[order(interesting_variants.stats$FDR, decreasing = FALSE), ]
    
# Generate the figures
plot_SV_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[1:9, ], 
                                    sample_info = sample_info, 
                                    singleton_column_name = c('Combined_unique_insertions'), 
                                    snp.mapping.info = variant.mapping.info,
                                    dir.output = dir.output, 
                                    output_name = paste('Singleton_Frequency_vs_Polymorphic_SVs_Set_1', sep = ''), 
                                    AFR_binary_genotypes = AFR_binary_genotypes, 
                                    AMR_binary_genotypes = AMR_binary_genotypes, 
                                    EUR_binary_genotypes = EUR_binary_genotypes, 
                                    EAS_binary_genotypes = EAS_binary_genotypes, 
                                    SAS_binary_genotypes = SAS_binary_genotypes,
                                    legend_position = 'bottomright') 

plot_SV_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[10:18, ], 
                                    sample_info = sample_info, 
                                    singleton_column_name = c('Combined_unique_insertions'), 
                                    snp.mapping.info = variant.mapping.info,
                                    dir.output = dir.output, 
                                    output_name = paste('Singleton_Frequency_vs_Polymorphic_SVs_Set_2', sep = ''), 
                                    AFR_binary_genotypes = AFR_binary_genotypes, 
                                    AMR_binary_genotypes = AMR_binary_genotypes, 
                                    EUR_binary_genotypes = EUR_binary_genotypes, 
                                    EAS_binary_genotypes = EAS_binary_genotypes, 
                                    SAS_binary_genotypes = SAS_binary_genotypes,
                                    legend_position = 'bottomright') 

plot_SV_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[19:21, ], 
                                    sample_info = sample_info, 
                                    singleton_column_name = c('Combined_unique_insertions'), 
                                    snp.mapping.info = variant.mapping.info,
                                    dir.output = dir.output, 
                                    output_name = paste('Singleton_Frequency_vs_Polymorphic_SVs_Set_3', sep = ''), 
                                    AFR_binary_genotypes = AFR_binary_genotypes, 
                                    AMR_binary_genotypes = AMR_binary_genotypes, 
                                    EUR_binary_genotypes = EUR_binary_genotypes, 
                                    EAS_binary_genotypes = EAS_binary_genotypes, 
                                    SAS_binary_genotypes = SAS_binary_genotypes,
                                    legend_position = 'bottomright') 










# PLOT GLOBAL SINGLETON FREQUENCY VS GENOTYPE FOR SNVS NEAR/OVERLAPPING KNOWN L1 TRANSPOSITION REGULATORS
    
# Define the snps of interest
interesting_variants <- c('rs1471205623', 'rs200381403', 'rs75237296', 'rs1333597286', 'rs1288384419')

# NOTES:
# rs200381403 - may overlap TRIM28 binding site
# rs75237296 - two COSMIC entries? one says FATHMM-MKL Score = 0.9919 (functionally significant/ detrimental)
# rs1333597286 - on COSMIC, FATHMM-MKL Score = 0.9918 (functionally significant/ detrimental)
    
# Extract the info for these snps
interesting_variants.stats <- sig.variants.greenlist[interesting_variants, c('SNP', 'Nearby_Genes', 'FDR', 'OR.R.')] # manually confirm that genes of interest are linked
    
    # Update the feature names associated with these snps
    interesting_variants.stats$Nearby_Genes <- c('MPHOSPH8', 'MPHOSPH8', 'PABPC1', 'PABPC1', 'RAD51B')
    
    # Add descriptor for the relationship between snps and genes
    interesting_variants.stats$Relationship <- c('upstream', 'upstream', 'in 3-UTR', 'in 3-UTR', 'intronic')
    
    # reorganize columns
    interesting_variants.stats <- interesting_variants.stats[, c('SNP', 'Nearby_Genes', 'Relationship', 'FDR', 'OR.R.')]
    
    # sort results by FDR
    interesting_variants.stats <- interesting_variants.stats[order(interesting_variants.stats$FDR, decreasing = FALSE), ]
    
# Generate the figures
plot_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[1:5, ], 
                                sample_info = sample_info, 
                                singleton_column_name = c('Combined_unique_insertions'), 
                                snp.mapping.info = variant.mapping.info,
                                dir.output = dir.output, 
                                output_name = paste('Singleton_Frequency_vs_SNVs_Near_Known_L1_Transposition_Regulators_Set_1', sep = ''), 
                                AFR_binary_genotypes = AFR_binary_genotypes, 
                                AMR_binary_genotypes = AMR_binary_genotypes, 
                                EUR_binary_genotypes = EUR_binary_genotypes, 
                                EAS_binary_genotypes = EAS_binary_genotypes, 
                                SAS_binary_genotypes = SAS_binary_genotypes,
                                legend_position = 'topright') 
      

    







# PLOT GLOBAL SINGLETON FREQUENCY VS GENOTYPE FOR SNVS NEAR/OVERLAPPING KNOWN L1 EXPRESSION REGULATORS
    
# Load variants near genes of interest
overlapping.expression.reg <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Enrichment_Results_SNVs/Greenlist_SNVs/Variants_overlapping_L1_Expression_Regulators.txt', header = TRUE, sep = "\t")

# Define the snps of interest
interesting_variants <- c('rs1350516110', 'rs201619112', 'rs71475866')
    
# Extract the info for these snps
interesting_variants.stats <- sig.variants.greenlist[interesting_variants, c('SNP', 'Nearby_Genes', 'FDR', 'OR.R.')] # manually confirm that genes of interest are linked
    
    # Update the feature names associated with these snps
    interesting_variants.stats$Nearby_Genes <- c('RHOT1', 'XPR1', 'PFKP')
    
    # Add descriptor for the relationship between snps and genes
    interesting_variants.stats$Relationship <- c('upstream', 'downstream', 'upstream')
    
    # reorganize columns
    interesting_variants.stats <- interesting_variants.stats[, c('SNP', 'Nearby_Genes', 'Relationship', 'FDR', 'OR.R.')]
    
    # sort results by FDR
    interesting_variants.stats <- interesting_variants.stats[order(interesting_variants.stats$FDR, decreasing = FALSE), ]
    
# Generate the figures
plot_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[1:3, ], 
                                sample_info = sample_info, 
                                singleton_column_name = c('Combined_unique_insertions'), 
                                snp.mapping.info = variant.mapping.info,
                                dir.output = dir.output, 
                                output_name = paste('Singleton_Frequency_vs_SNVs_Near_Known_L1_Expression_Regulators_Set_1', sep = ''), 
                                AFR_binary_genotypes = AFR_binary_genotypes, 
                                AMR_binary_genotypes = AMR_binary_genotypes, 
                                EUR_binary_genotypes = EUR_binary_genotypes, 
                                EAS_binary_genotypes = EAS_binary_genotypes, 
                                SAS_binary_genotypes = SAS_binary_genotypes,
                                legend_position = 'topright') 
      

    







# PLOT GLOBAL SINGLETON FREQUENCY VS GENOTYPE FOR SNVS NEAR/OVERLAPPING METHYLTRANSFERASES
    
# Load variants near genes of interest
overlapping.methyltransferase <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Enrichment_Results_SNVs/Greenlist_SNVs/Variants_overlapping_Histone_Methyltransferases.txt', header = TRUE, sep = "\t")

# Define the snps of interest
interesting_variants <- c('rs369852550', 'rs375298724')
    
# Extract the info for these snps
interesting_variants.stats <- sig.variants.greenlist[interesting_variants, c('SNP', 'Nearby_Genes', 'FDR', 'OR.R.')] # manually confirm that genes of interest are linked
    
    # Update the feature names associated with these snps
    interesting_variants.stats$Nearby_Genes <- c('EEF2KMT', 'PRDM7')
    
    # Add descriptor for the relationship between snps and genes
    interesting_variants.stats$Relationship <- c('upstream', 'upstream')
    
    # reorganize columns
    interesting_variants.stats <- interesting_variants.stats[, c('SNP', 'Nearby_Genes', 'Relationship', 'FDR', 'OR.R.')]
    
    # sort results by FDR
    interesting_variants.stats <- interesting_variants.stats[order(interesting_variants.stats$FDR, decreasing = FALSE), ]
    
# Generate the figures
plot_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[1:2, ], 
                                sample_info = sample_info, 
                                singleton_column_name = c('Combined_unique_insertions'), 
                                snp.mapping.info = variant.mapping.info,
                                dir.output = dir.output, 
                                output_name = paste('Singleton_Frequency_vs_SNVs_Near_Methyltransferases_Set_1', sep = ''), 
                                AFR_binary_genotypes = AFR_binary_genotypes, 
                                AMR_binary_genotypes = AMR_binary_genotypes, 
                                EUR_binary_genotypes = EUR_binary_genotypes, 
                                EAS_binary_genotypes = EAS_binary_genotypes, 
                                SAS_binary_genotypes = SAS_binary_genotypes,
                                legend_position = 'topright') 
      

    







# PLOT GLOBAL SINGLETON FREQUENCY VS GENOTYPE FOR SNVS NEAR/OVERLAPPING RNA MODIFICATION GENES
    
# Load variants near genes of interest
overlapping.RNAmod <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Enrichment_Results_SNVs/Greenlist_SNVs/Variants_overlapping_RNA_modification_genes.txt', header = TRUE, sep = "\t")

# Define the snps of interest
interesting_variants <- c('rs71475866', 'rs865777935', 'rs144587251')
    
# Extract the info for these snps
interesting_variants.stats <- sig.variants.greenlist[interesting_variants, c('SNP', 'Nearby_Genes', 'FDR', 'OR.R.')] # manually confirm that genes of interest are linked
    
    # Update the feature names associated with these snps
    interesting_variants.stats$Nearby_Genes <- c('ADARB2', 'A1CF', 'METTL14')
    
    # Add descriptor for the relationship between snps and genes
    interesting_variants.stats$Relationship <- c('upstream', 'downstream', 'upstream')
    
    # reorganize columns
    interesting_variants.stats <- interesting_variants.stats[, c('SNP', 'Nearby_Genes', 'Relationship', 'FDR', 'OR.R.')]
    
    # sort results by FDR
    interesting_variants.stats <- interesting_variants.stats[order(interesting_variants.stats$FDR, decreasing = FALSE), ]
    
# Generate the figures
plot_genotype_vs_Singleton_Freq(snp_feature_df = interesting_variants.stats[1:3, ], 
                                sample_info = sample_info, 
                                singleton_column_name = c('Combined_unique_insertions'), 
                                snp.mapping.info = variant.mapping.info,
                                dir.output = dir.output, 
                                output_name = paste('Singleton_Frequency_vs_SNVs_Near_RNA_Modification_Genes_Set_1', sep = ''), 
                                AFR_binary_genotypes = AFR_binary_genotypes, 
                                AMR_binary_genotypes = AMR_binary_genotypes, 
                                EUR_binary_genotypes = EUR_binary_genotypes, 
                                EAS_binary_genotypes = EAS_binary_genotypes, 
                                SAS_binary_genotypes = SAS_binary_genotypes,
                                legend_position = 'topright') 





# Clean the environment
rm(list=ls())       
        
