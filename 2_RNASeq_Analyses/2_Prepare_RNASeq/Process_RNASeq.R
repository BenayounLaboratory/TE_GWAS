# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Process_RNASeq_functions.R")

# Load libraries
library(edgeR)
library(genefilter)
library(DESeq2) # For VST
library(BiocParallel) # Used by DESeq2 for parallelization
  register(MulticoreParam(6)) # Use six cores
library(limma) # for batch effect removal
  
# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/'

    
    
# Section 1: FILTER LOW EXPRESSION GENES ------------------------------------------------------------------------------------------------------------------------


  
# Load sample meta-data
GEUV_sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Comment.ENA_RUN.

# Define directories and lists for RNAseq count tables.
count_Table.dir.YRI <- c("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/Quantifications_from_Bravo_et_al_2024/TETranscripts_Counts_YRI/")
count_Table.dir.EUR <- c("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/Quantifications_from_Bravo_et_al_2024/TETranscripts_Counts_EUR/")

cntTable_complete.list.YRI <- list.files(count_Table.dir.YRI, "\\.cntTable$", recursive=FALSE, full.names=TRUE)
cntTable_complete.list.EUR <- list.files(count_Table.dir.EUR, "\\.cntTable$", recursive=FALSE, full.names=TRUE)

# Run function to combine the individual counts tables.
gene_TE_counts_raw_YRI <- aggregate_counts(cntTable_complete.list = cntTable_complete.list.YRI)
gene_TE_counts_raw_EUR <- aggregate_counts(cntTable_complete.list = cntTable_complete.list.EUR)

# Shorten sample column names 
colnames(gene_TE_counts_raw_YRI) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw_YRI))
colnames(gene_TE_counts_raw_EUR) <- sub("Aligned.sortedByCoord.out.bam*", "", colnames(gene_TE_counts_raw_EUR))

colnames(gene_TE_counts_raw_YRI) <- sub("*fastp_", "", colnames(gene_TE_counts_raw_YRI))
colnames(gene_TE_counts_raw_EUR) <- sub("*fastp_", "", colnames(gene_TE_counts_raw_EUR))

# change column names to match the names used in the genotyping files.
colnames(gene_TE_counts_raw_YRI) <- GEUV_sample_info[colnames(gene_TE_counts_raw_YRI), 1]
colnames(gene_TE_counts_raw_EUR) <- GEUV_sample_info[colnames(gene_TE_counts_raw_EUR), 1]

# Run function to remove gene version info, combine same genes, and filter lowly expressed genes. *IMPORTANT: PSEUDOAUTOSOMAL GENES (PAR_Y) ARE COMBINED
gene_TE_counts.YRI <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_YRI, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 0.50) 
gene_TE_counts.EUR <- cleanup_and_filter_counts(count_data = gene_TE_counts_raw_EUR, filter = TRUE, min_counts_per_sample = 10, fraction_of_samples = 0.50) 

    # Find genes in common between the two populations (for generating consensus networks later)
    common_genes <- rownames(gene_TE_counts.YRI[rownames(gene_TE_counts.YRI) %in% rownames(gene_TE_counts.EUR), ])
    
    # Subset common genes
    gene_TE_counts.YRI <- gene_TE_counts.YRI[common_genes, ]
    gene_TE_counts.EUR <- gene_TE_counts.EUR[common_genes, ]

    # Save counts files.
    write.table(gene_TE_counts.YRI, file = paste(dir.output, "All_counts_YRI_86_filtered", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    write.table(gene_TE_counts.EUR, file = paste(dir.output, "All_counts_EUR_358_filtered", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = FALSE)
    


# Section 2: BATCH EFFECT REMOVAL -------------------------------------

      
      
# DEFINE KNOWN COVARIATES
      
# Re-load sample meta data
GEUV_sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Source.Name

# Load genotype PCA covariates
Principal_components <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/1_Population_Structure_Analysis_SNV_and_SV/Population_Structure_Analysis/COVARIATE_Population_Structure_PCA_SNVs_and_SVs.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load the L1/Alu insertion/deletion covariate table
# L1_Alu_insertions_deletions <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_eQTL_Analysis/1_Prepare_Genotype_Data/Structural_Variant_Counts/L1_Alu_Insertions_Deletions_Per_Sample.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
    
# Extract known covariates from the meta data
covariates_table <- GEUV_sample_info[, c('Source.Name','Characteristics.sex.', 'Characteristics.ancestry.category.', 'Performer')]

# Update covariate names
colnames(covariates_table) <- c('id', 'sex', 'ancestry', 'lab') 

# Fill in SNV-SV genotype PCs 1-2 (note: these were calculated separately in each superpopulation)
covariates_table$PC1 <- Principal_components[rownames(covariates_table), 'PC1']
covariates_table$PC2 <- Principal_components[rownames(covariates_table), 'PC2']

# Add the net L1/Alu copy number (ie Alu/L1 insertions - Alu/L1 deletions)
# covariates_table$Net_L1_Alu_Copies <- L1_Alu_insertions_deletions[rownames(covariates_table), 'Net_Both']

# Split covariates into YRI and EUR
covariates_table.YRI <- covariates_table[which(covariates_table$ancestry == 'Yoruba'), ]
covariates_table.EUR <- covariates_table[which(covariates_table$ancestry != 'Yoruba'), ]

# Scale and center the net copy number in each population
# covariates_table.YRI$Net_L1_Alu_Copies <- scale(covariates_table.YRI$Net_L1_Alu_Copies, center = TRUE, scale = TRUE)
# covariates_table.EUR$Net_L1_Alu_Copies <- scale(covariates_table.EUR$Net_L1_Alu_Copies, center = TRUE, scale = TRUE)

# Order expression data samples like the covariates table
gene_TE_counts.YRI <- gene_TE_counts.YRI[, rownames(covariates_table.YRI)]
gene_TE_counts.EUR <- gene_TE_counts.EUR[, rownames(covariates_table.EUR)]


    

    
    
    
# VST TRANSFORMATION OF THE COUNT DATA

# Collect covariates in a dataframe
SampleInfo.YRI <- data.frame(row.names = rownames(covariates_table.YRI),
                             sex = covariates_table.YRI$sex,
                             lab = covariates_table.YRI$lab,
                             PC1 = covariates_table.YRI$PC1,
                             PC2 = covariates_table.YRI$PC2)
                             
SampleInfo.EUR <- data.frame(row.names = rownames(covariates_table.EUR),
                             sex = covariates_table.EUR$sex,
                             ancestry = covariates_table.EUR$ancestry,
                             lab = covariates_table.EUR$lab,
                             PC1 = covariates_table.EUR$PC1,
                             PC2 = covariates_table.EUR$PC2) 
  
# Create DESeq2 object 
dds.YRI <- DESeqDataSetFromMatrix(countData = gene_TE_counts.YRI,
                                  colData = SampleInfo.YRI,
                                  design = ~ sex + lab + PC1 + PC2)

dds.EUR <- DESeqDataSetFromMatrix(countData = gene_TE_counts.EUR,
                                  colData = SampleInfo.EUR,
                                  design = ~ sex + lab + ancestry + PC1 + PC2)
  
# Run DESeq2
dds.YRI <- DESeq(dds.YRI, parallel = TRUE)
dds.EUR <- DESeq(dds.EUR, parallel = TRUE)

# Get VST data (includes size factor normalization)
VST.YRI <- assay(varianceStabilizingTransformation(dds.YRI, blind = FALSE))
VST.EUR <- assay(varianceStabilizingTransformation(dds.EUR, blind = FALSE))

    # Save VST data
    write.table(VST.YRI, file = paste(dir.output, 'All_counts_YRI_86_filtered_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
    write.table(VST.EUR, file = paste(dir.output, 'All_counts_EUR_358_filtered_VST', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
        
    
    
    
    
    
    
# OBTAIN EBV EXPRESSION AS AN ADDITIONAL COVARIATE
    
# Add EBV VST expression to the covariates table    
covariates_table.YRI$EBV_expr_VST <- VST.YRI["EBV_gene", rownames(covariates_table.YRI)]
covariates_table.EUR$EBV_expr_VST <- VST.EUR["EBV_gene", rownames(covariates_table.EUR)]
    
# Scale and center EBV expression
covariates_table.YRI$EBV_expr_VST <- scale(covariates_table.YRI$EBV_expr_VST, center = TRUE, scale = TRUE)
covariates_table.EUR$EBV_expr_VST <- scale(covariates_table.EUR$EBV_expr_VST, center = TRUE, scale = TRUE)
    
# Remove redundant variables
covariates_table.YRI <- covariates_table.YRI[ , -c(3)] # NOTE: DOUBLE CHECK THE COLUMNS REMOVED. REMOVE ANCESTRY.
covariates_table.EUR <- covariates_table.EUR[, ] # NOTE: DOUBLE CHECK COLUMNS REMOVED.

# Save covariates tables
write.table(covariates_table.YRI, file = paste(dir.output, "COVARIATES_YRI", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)
write.table(covariates_table.EUR, file = paste(dir.output, "COVARIATES_EUR", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)

    

    
    
    

# BATCH EFFECT REMOVAL WITH LIMMA
    
# Make model matrices
full.model.YRI = model.matrix(~ sex + lab + PC1 + PC2 + EBV_expr_VST, data = covariates_table.YRI)
full.model.EUR = model.matrix(~ sex + lab + ancestry + PC1 + PC2 + EBV_expr_VST, data = covariates_table.EUR)
           
# From the VST data, regress out lab + ancestry + PCs + sex + EBV
batch_corr_YRI <- removeBatchEffect(x = VST.YRI,
                                    batch = NULL,
                                    covariates = full.model.YRI[, c(2:11)], 
                                    design = full.model.YRI[, c(1), drop = FALSE])

batch_corr_EUR <- removeBatchEffect(x = VST.EUR,
                                    batch = NULL,
                                    covariates = full.model.EUR[, c(2:14)], 
                                    design = full.model.EUR[, c(1), drop = FALSE])

# Remove EBV from the expression matrix, since variation has been removed

    # Find indices for EBV
    EBV.index.YRI <- which(rownames(batch_corr_YRI) %in% 'EBV_gene')
    EBV.index.EUR <- which(rownames(batch_corr_EUR) %in% 'EBV_gene')
    
    # Remove EBV
    batch_corr_YRI <- batch_corr_YRI[-c(EBV.index.YRI), ]
    batch_corr_EUR <- batch_corr_EUR[-c(EBV.index.EUR), ]

# Save batch corrected data
write.table(batch_corr_YRI, file = paste(dir.output, 'All_counts_YRI_86_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_EBV', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)
write.table(batch_corr_EUR, file = paste(dir.output, 'All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_EBV', '.txt', sep =""), sep = "\t", row.names = T, col.names = NA, quote = F)

            
            
# SESSION INFO ---------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Process_RNASeq.txt", sep =""))
sessionInfo()
sink()    
    

# Clean the environment
rm(list=ls())

