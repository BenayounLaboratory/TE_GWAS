# Set strings as factors
options(stringsAsFactors = F)


# PREPARE COVARIATES AND PHENOTYPE DATA FOR PLINK ---------------------------
# The files prepared here will be used to run the GWAS


# Define output directories
dir.output.AFR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/'
dir.output.AMR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_Genotypes/'
dir.output.EAS <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_Genotypes/'
dir.output.EUR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes/'
dir.output.SAS <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_Genotypes/'

# Load sample metadata with GWAS case/control status
sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load SNV+SV PCA covariates
Pop_structure.PCA <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Population_Structure_Analysis/SNV_AND_SV_Population_Structure_PCA.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Update sample order to match sample metadata
    Pop_structure.PCA <- Pop_structure.PCA[rownames(sample_info), ]

# Append SNV+SV PCs 1-4 to the sample data
sample_info$PC1 <- Pop_structure.PCA$PC1
sample_info$PC2 <- Pop_structure.PCA$PC2
sample_info$PC3 <- Pop_structure.PCA$PC3
sample_info$PC4 <- Pop_structure.PCA$PC4

# Subset samples by superpopulation
samples.AFR <- sample_info[sample_info$Super_Population == 'AFR', ]
samples.AMR <- sample_info[sample_info$Super_Population == 'AMR', ]
samples.EAS <- sample_info[sample_info$Super_Population == 'EAS', ]
samples.EUR <- sample_info[sample_info$Super_Population == 'EUR', ]
samples.SAS <- sample_info[sample_info$Super_Population == 'SAS', ]

# Make model matrices with covariates and phenotypes of interest 
# Note 1: The GWAS will primarily focus on the combined Alu+L1 singleton groups. I will explore associations using the combined target site duplication (TSD)- containing Alu+L1 group
# Note 2: Since I'm correcting for population structure though both SNV and SV PCs, I'm not going to include the population label, since it would be colinear with the PCs and PLINK gives a collinarity warning
# Note 3: Using model matrices is necessary since categorical variables (like sex) need to be in numerical format for Plink
samples.AFR <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + Sex + Combined_Plink_Binary + Combined_TSD_Plink_Binary, data = samples.AFR) 
samples.AMR <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + Sex + Combined_Plink_Binary + Combined_TSD_Plink_Binary, data = samples.AMR) 
samples.EAS <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + Sex + Combined_Plink_Binary + Combined_TSD_Plink_Binary, data = samples.EAS) 
samples.EUR <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + Sex + Combined_Plink_Binary + Combined_TSD_Plink_Binary, data = samples.EUR) 
samples.SAS <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + Sex + Combined_Plink_Binary + Combined_TSD_Plink_Binary, data = samples.SAS) 

# Remove intercept column
samples.AFR <- samples.AFR[, -c(1)] 
samples.AMR <- samples.AMR[, -c(1)] 
samples.EAS <- samples.EAS[, -c(1)] 
samples.EUR <- samples.EUR[, -c(1)] 
samples.SAS <- samples.SAS[, -c(1)]

# Convert matrices back to dataframes
samples.AFR <- as.data.frame(samples.AFR)
samples.AMR <- as.data.frame(samples.AMR)
samples.EAS <- as.data.frame(samples.EAS)
samples.EUR <- as.data.frame(samples.EUR)
samples.SAS <- as.data.frame(samples.SAS)

# Update male/female values to reflect Plink encodings
# 1 or M designates males in Plink
# 2 or F designates females in Plink
samples.AFR[samples.AFR$Sexmale == 1, 'Sexmale'] <- 1 
samples.AFR[samples.AFR$Sexmale == 0, 'Sexmale'] <- 2 

samples.AMR[samples.AMR$Sexmale == 1, 'Sexmale'] <- 1 
samples.AMR[samples.AMR$Sexmale == 0, 'Sexmale'] <- 2 
      
samples.EAS[samples.EAS$Sexmale == 1, 'Sexmale'] <- 1 
samples.EAS[samples.EAS$Sexmale == 0, 'Sexmale'] <- 2 
      
samples.EUR[samples.EUR$Sexmale == 1, 'Sexmale'] <- 1 
samples.EUR[samples.EUR$Sexmale == 0, 'Sexmale'] <- 2 

samples.SAS[samples.SAS$Sexmale == 1, 'Sexmale'] <- 1
samples.SAS[samples.SAS$Sexmale == 0, 'Sexmale'] <- 2 
      
# Generate a file with sample names to keep in each analysis

    # Define samples names
    keep.AFR <- cbind(rownames(samples.AFR), rownames(samples.AFR))
    keep.AMR <- cbind(rownames(samples.AMR), rownames(samples.AMR))
    keep.EAS <- cbind(rownames(samples.EAS), rownames(samples.EAS))
    keep.EUR <- cbind(rownames(samples.EUR), rownames(samples.EUR))
    keep.SAS <- cbind(rownames(samples.SAS), rownames(samples.SAS))
      
    # Save files
    write.table(keep.AFR, file = paste(dir.output.AFR, "Keep_AFR_660", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(keep.AMR, file = paste(dir.output.AMR, "Keep_AMR_347", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(keep.EAS, file = paste(dir.output.EAS, "Keep_EAS_504", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(keep.EUR, file = paste(dir.output.EUR, "Keep_EUR_503", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(keep.SAS, file = paste(dir.output.SAS, "Keep_SAS_489", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)

# Generate a file specifying sample sex info

    # Define sample sexes
    sex.AFR <- cbind(rownames(samples.AFR), rownames(samples.AFR), samples.AFR$Sexmale)
    sex.AMR <- cbind(rownames(samples.AMR), rownames(samples.AMR), samples.AMR$Sexmale)
    sex.EAS <- cbind(rownames(samples.EAS), rownames(samples.EAS), samples.EAS$Sexmale)
    sex.EUR <- cbind(rownames(samples.EUR), rownames(samples.EUR), samples.EUR$Sexmale)
    sex.SAS <- cbind(rownames(samples.SAS), rownames(samples.SAS), samples.SAS$Sexmale)
      
    # Save files
    write.table(sex.AFR, file = paste(dir.output.AFR, "Sex_AFR_660", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(sex.AMR, file = paste(dir.output.AMR, "Sex_AMR_347", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(sex.EAS, file = paste(dir.output.EAS, "Sex_EAS_504", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(sex.EUR, file = paste(dir.output.EUR, "Sex_EUR_503", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(sex.SAS, file = paste(dir.output.SAS, "Sex_SAS_489", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)

# Generate a file specifying covariates  
# NOTE: If VIF correction is needed, you can use the "multicol" function in "fuzzySim" package

    # Define covariates (sex will be specified in a separate file, so only SNV+SV PCs here)
    covariates.AFR <- cbind(rownames(samples.AFR), rownames(samples.AFR), samples.AFR[, 1:4])
    covariates.AMR <- cbind(rownames(samples.AMR), rownames(samples.AMR), samples.AMR[, 1:4])
    covariates.EAS <- cbind(rownames(samples.EAS), rownames(samples.EAS), samples.EAS[, 1:4])
    covariates.EUR <- cbind(rownames(samples.EUR), rownames(samples.EUR), samples.EUR[, 1:4]) 
    covariates.SAS <- cbind(rownames(samples.SAS), rownames(samples.SAS), samples.SAS[, 1:4])

    # Save files
    write.table(covariates.AFR, file = paste(dir.output.AFR, "Covariates_AFR_660", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(covariates.AMR, file = paste(dir.output.AMR, "Covariates_AMR_347", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(covariates.EAS, file = paste(dir.output.EAS, "Covariates_EAS_504", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(covariates.EUR, file = paste(dir.output.EUR, "Covariates_EUR_503", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(covariates.SAS, file = paste(dir.output.SAS, "Covariates_SAS_489", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)

# Generate a file specifying phenotypes (presence or absence of global L1/Alu singletons)

    # Define phenotypes
    phenotype.AFR <- cbind(rownames(samples.AFR), rownames(samples.AFR), samples.AFR$Combined_Plink_Binary)
    phenotype.AMR <- cbind(rownames(samples.AMR), rownames(samples.AMR), samples.AMR$Combined_Plink_Binary)
    phenotype.EAS <- cbind(rownames(samples.EAS), rownames(samples.EAS), samples.EAS$Combined_Plink_Binary)
    phenotype.EUR <- cbind(rownames(samples.EUR), rownames(samples.EUR), samples.EUR$Combined_Plink_Binary)
    phenotype.SAS <- cbind(rownames(samples.SAS), rownames(samples.SAS), samples.SAS$Combined_Plink_Binary)
    
    # Save files
    write.table(phenotype.AFR, file = paste(dir.output.AFR, "Phenotype_Combined_L1_Alu_Singleton_AFR_660", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(phenotype.AMR, file = paste(dir.output.AMR, "Phenotype_Combined_L1_Alu_Singleton_AMR_347", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(phenotype.EAS, file = paste(dir.output.EAS, "Phenotype_Combined_L1_Alu_Singleton_EAS_504", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(phenotype.EUR, file = paste(dir.output.EUR, "Phenotype_Combined_L1_Alu_Singleton_EUR_503", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(phenotype.SAS, file = paste(dir.output.SAS, "Phenotype_Combined_L1_Alu_Singleton_SAS_489", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    
# Generate a file specifying secondary phenotypes (presence or absence of global L1/Alu singletons WITH TARGET SITE DUPLICATIONS)

    # Define phenotypes
    TSD.phenotype.AFR <- cbind(rownames(samples.AFR), rownames(samples.AFR), samples.AFR$Combined_TSD_Plink_Binary)
    TSD.phenotype.AMR <- cbind(rownames(samples.AMR), rownames(samples.AMR), samples.AMR$Combined_TSD_Plink_Binary)
    TSD.phenotype.EAS <- cbind(rownames(samples.EAS), rownames(samples.EAS), samples.EAS$Combined_TSD_Plink_Binary)
    TSD.phenotype.EUR <- cbind(rownames(samples.EUR), rownames(samples.EUR), samples.EUR$Combined_TSD_Plink_Binary)
    TSD.phenotype.SAS <- cbind(rownames(samples.SAS), rownames(samples.SAS), samples.SAS$Combined_TSD_Plink_Binary)
    
    # Save files
    write.table(TSD.phenotype.AFR, file = paste(dir.output.AFR, "Phenotype_TSD_Combined_L1_Alu_Singleton_AFR_660", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(TSD.phenotype.AMR, file = paste(dir.output.AMR, "Phenotype_TSD_Combined_L1_Alu_Singleton_AMR_347", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(TSD.phenotype.EAS, file = paste(dir.output.EAS, "Phenotype_TSD_Combined_L1_Alu_Singleton_EAS_504", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(TSD.phenotype.EUR, file = paste(dir.output.EUR, "Phenotype_TSD_Combined_L1_Alu_Singleton_EUR_503", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)
    write.table(TSD.phenotype.SAS, file = paste(dir.output.SAS, "Phenotype_TSD_Combined_L1_Alu_Singleton_SAS_489", ".txt", sep=""), row.names = FALSE, col.names = FALSE, sep = "\t", na = "NA", quote = F)

  
      
      
      
# Remove unneeded variables
rm(list=ls()) 
      
      
