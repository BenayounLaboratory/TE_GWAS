# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/WGCNA_Functions.R") 

# Load libraries
library(WGCNA) # for constructing the gene network
library(ggplot2) # for making plots
library(cowplot) # for ggplot image grids?
library(BEDMatrix)
library(poolr) # to combine pvalues

# Define output directory
dir.trait.cor <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Trait_Correlations_Consensus/'



# LOAD DATA --------------------------------------------


# DEFINE INPUT PARAMETERS

# Define the number of datasets
nSets = 2

# Load expression data.
VST.expression.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_EBV.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
VST.expression.YRI <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/All_counts_YRI_86_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_EBV.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Transpose the data to match WGCNA input requirements
    VST.expression.EUR <- t(VST.expression.EUR)
    VST.expression.YRI <- t(VST.expression.YRI)
    
    # Define the number of samples
    my.nSamples.EUR = nrow(VST.expression.EUR)
    my.nSamples.YRI = nrow(VST.expression.YRI)
    
# Put the expression data into a multiset for consensus network generation. 
my.multiExpr = vector(mode = "list", length = nSets)

    my.multiExpr[[1]] = list(data = as.data.frame(VST.expression.EUR))
    names(my.multiExpr[[1]]$data) = colnames(VST.expression.EUR)
    rownames(my.multiExpr[[1]]$data) = rownames(VST.expression.EUR)
    
    my.multiExpr[[2]] = list(data = as.data.frame(VST.expression.YRI))
    names(my.multiExpr[[2]]$data) = colnames(VST.expression.YRI)
    rownames(my.multiExpr[[2]]$data) = rownames(VST.expression.YRI)
    
    # Check that the data has the correct format for many functions operating on multiple sets:
    checkSets(my.multiExpr)
    
# Load genotype data
    
    # Specify the genotype BED files
    plink_file_path.YRI <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/Final_SNVs_SVs_AFR.bed'
    plink_file_path.EUR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes/Final_SNVs_SVs_EUR.bed'
    
    # Make R object to stream snp genotypes into memory (rows are samples and columns are snps). # NOTE: If used, MT snps will need to be assigned a snpid, since they dont have rsid.
    binary_genotypes.YRI <- BEDMatrix(path = plink_file_path.YRI, simple_names = T)
    binary_genotypes.EUR <- BEDMatrix(path = plink_file_path.EUR, simple_names = T)

# Load the results of network analysis 
load(file = "/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Network_Construction_Consensus/TETranscripts_GEUVADIS_CONSENSUS_MEs_and_moduleColors.RData")

    # Define module color names
    MEColorNames <- colnames(my.consMEs[[1]][["data"]])
  
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("EUR LCLs", "YRI LCLs")

    


    
# CORRELATION TO CONTROL/CASE STATUS (L1 vs Alu vs ALL) --------------------------------------------
     
    
# Load case/control status
singleton.counts <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Use the binary case/control columns to distinguish controls (1) from cases (2)
    # Subtract 1 from those columns so that controls are 0 and cases are 1
    singleton.counts$L1_Plink_Binary <- singleton.counts$L1_Plink_Binary - 1
    singleton.counts$Alu_Plink_Binary <- singleton.counts$Alu_Plink_Binary - 1
    singleton.counts$Combined_Plink_Binary <- singleton.counts$Combined_Plink_Binary - 1

    # Only keep samples with transcriptomes
    singleton.counts.EUR <- singleton.counts[rownames(VST.expression.EUR), ]
    singleton.counts.YRI <- singleton.counts[rownames(VST.expression.YRI), ]

# Define the phenotypes (case/control status for L1, Alu, or Either singleton status)
my.phenotype.EUR <- singleton.counts.EUR[, c('L1_Plink_Binary', 'Alu_Plink_Binary', 'Combined_Plink_Binary')]
my.phenotype.YRI <- singleton.counts.YRI[, c('L1_Plink_Binary', 'Alu_Plink_Binary', 'Combined_Plink_Binary')]

# Put the phenotype data into a multiset 
my.Traits = vector(mode = "list", length = nSets)

    my.Traits[[1]] = list(data = as.data.frame(my.phenotype.EUR))
    names(my.Traits[[1]]$data) = colnames(my.phenotype.EUR)
    rownames(my.Traits[[1]]$data) = rownames(my.phenotype.EUR)
    
    my.Traits[[2]] = list(data = as.data.frame(my.phenotype.YRI))
    names(my.Traits[[2]]$data) = colnames(my.phenotype.YRI)
    rownames(my.Traits[[2]]$data) = rownames(my.phenotype.YRI)

# Set up variables to contain the module-trait correlations/pvalues
my.Trait.Cor = list()
my.Trait.Pvalue = list()

# Calculate the correlations/pvalues
for (set in 1:nSets){
  
    my.Trait.Cor[[set]] = cor(my.consMEs[[set]]$data, my.Traits[[set]]$data, use = "p")
    
    my.Trait.Pvalue[[set]] = corPvalueFisher(my.Trait.Cor[[set]], nrow(my.consMEs[[set]][["data"]]))
    
}


# MAKE PLOT (SET1)

    # Define the set number
    set = 1

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_L1_Alu_Combined_Case_Control_Status_EUR.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Case/Control Status in", setLabels[set]))
    
    dev.off()

    
# MAKE PLOT (SET2)

    # Define the set number
    set = 2

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_L1_Alu_Combined_Case_Control_Status_YRI.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Case/Control Status in", setLabels[set]))
    
    dev.off()





# COMBINE RESULTS
    
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))
consensusPvalue = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))

# Find consensus negative correlations
negative = my.Trait.Cor[[1]] < 0 & my.Trait.Cor[[2]] < 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = negative, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = negative, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

# Find consensus positive correlations
positive = my.Trait.Cor[[1]] > 0 & my.Trait.Cor[[2]] > 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = positive, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = positive, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

    
# MAKE PLOT (META-ANALYSIS)

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(consensusCor, 2), "\n(",
                        signif(consensusPvalue, 1), ")", sep = "")

    dim(textMatrix) = dim(consensusCor)

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_L1_Alu_Combined_Case_Control_Status_META.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = consensusCor,
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Case/Control Status (EUR & YRI)"))
    
    dev.off()


    
    
                 
# CORRELATION TO CONTROL/CASE STATUS (COMBINED ONLY) --------------------------------------------
     
    
# Load case/control status
singleton.counts <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Use the binary case/control columns to distinguish controls (1) from cases (2)
    # Subtract 1 from those columns so that controls are 0 and cases are 1
    singleton.counts$Combined_Plink_Binary <- singleton.counts$Combined_Plink_Binary - 1

    # Only keep samples with transcriptomes
    singleton.counts.EUR <- singleton.counts[rownames(VST.expression.EUR), ]
    singleton.counts.YRI <- singleton.counts[rownames(VST.expression.YRI), ]

# Define the phenotypes (case/control status for L1, Alu, or Either singleton status)
my.phenotype.EUR <- singleton.counts.EUR[, c('Combined_Plink_Binary'), drop = FALSE]
my.phenotype.YRI <- singleton.counts.YRI[, c('Combined_Plink_Binary'), drop = FALSE]

# Put the phenotype data into a multiset 
my.Traits = vector(mode = "list", length = nSets)

    my.Traits[[1]] = list(data = as.data.frame(my.phenotype.EUR))
    names(my.Traits[[1]]$data) = colnames(my.phenotype.EUR)
    rownames(my.Traits[[1]]$data) = rownames(my.phenotype.EUR)
    
    my.Traits[[2]] = list(data = as.data.frame(my.phenotype.YRI))
    names(my.Traits[[2]]$data) = colnames(my.phenotype.YRI)
    rownames(my.Traits[[2]]$data) = rownames(my.phenotype.YRI)

# Set up variables to contain the module-trait correlations/pvalues
my.Trait.Cor = list()
my.Trait.Pvalue = list()

# Calculate the correlations/pvalues
for (set in 1:nSets){
  
    my.Trait.Cor[[set]] = cor(my.consMEs[[set]]$data, my.Traits[[set]]$data, use = "p")
    
    my.Trait.Pvalue[[set]] = corPvalueFisher(my.Trait.Cor[[set]], nrow(my.consMEs[[set]][["data"]]))
    
}


# MAKE PLOT (FOR EUR CONSENSUS NETWORK)

    # Define the set number
    set = 1

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Combined_Case_Control_Status_EUR.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Case/Control Status in", setLabels[set]))
    
    dev.off()

    
# MAKE PLOT (FOR YRI CONSENSUS NETWORK)

    # Define the set number
    set = 2

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Case_Combined_Control_Status_YRI.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Case/Control Status in", setLabels[set]))
    
    dev.off()





# COMBINE RESULTS
    
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))
consensusPvalue = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))

# Find consensus negative correlations
negative = my.Trait.Cor[[1]] < 0 & my.Trait.Cor[[2]] < 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = negative, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = negative, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

# Find consensus positive correlations
positive = my.Trait.Cor[[1]] > 0 & my.Trait.Cor[[2]] > 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = positive, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = positive, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

    
# MAKE PLOT (META-ANALYSIS)

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(consensusCor, 2), "\n(",
                        signif(consensusPvalue, 1), ")", sep = "")

    dim(textMatrix) = dim(consensusCor)

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Combined_Case_Control_Status_META.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = consensusCor,
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Case/Control Status (EUR & YRI)"))
    
    dev.off()


    
    
                 
# CORRELATION TO VARIANTS NEAR/OVERLAPPING L1 TRANSPOSITION REGULATORS --------------------------------------------
     
    

# Define the variants of interest
my.variants <- c('rs75237296', 'rs1288384419', 'rs1471205623') 

# Define plot labels for variants
my.variant.labels <- c('PABPC1', 'RAD51B', 'MPHOSPH8')

# Extract genotypes and define them as the phenotype
my.phenotype.EUR <- as.data.frame(binary_genotypes.EUR[rownames(my.consMEs[[1]][["data"]]), my.variants])
my.phenotype.YRI <- as.data.frame(binary_genotypes.YRI[rownames(my.consMEs[[2]][["data"]]), my.variants])

# Put the phenotype data into a multiset 
my.Traits = vector(mode = "list", length = nSets)

    my.Traits[[1]] = list(data = as.data.frame(my.phenotype.EUR))
    names(my.Traits[[1]]$data) = colnames(my.phenotype.EUR)
    rownames(my.Traits[[1]]$data) = rownames(my.phenotype.EUR)
    
    my.Traits[[2]] = list(data = as.data.frame(my.phenotype.YRI))
    names(my.Traits[[2]]$data) = colnames(my.phenotype.YRI)
    rownames(my.Traits[[2]]$data) = rownames(my.phenotype.YRI)

# Set up variables to contain the module-trait correlations/pvalues
my.Trait.Cor = list()
my.Trait.Pvalue = list()

# Calculate the correlations/pvalues
for (set in 1:nSets){
  
    my.Trait.Cor[[set]] = cor(my.consMEs[[set]]$data, my.Traits[[set]]$data, use = "p")
    
    my.Trait.Pvalue[[set]] = corPvalueFisher(my.Trait.Cor[[set]], nrow(my.consMEs[[set]][["data"]]))
    
}


# MAKE PLOT (SET1)

    # Define the set number
    set = 1

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Transposition_Regulator_SNVs_EUR.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with L1 Regulators SNVs in", setLabels[set]))
    
    dev.off()

    
# MAKE PLOT (SET2)

    # Define the set number
    set = 2

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Transposition_Regulator_SNVs_YRI.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with L1 Regulators SNVs in", setLabels[set]))
    
    dev.off()





# COMBINE RESULTS
    
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))
consensusPvalue = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))

# Find consensus negative correlations
negative = my.Trait.Cor[[1]] < 0 & my.Trait.Cor[[2]] < 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = negative, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = negative, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

# Find consensus positive correlations
positive = my.Trait.Cor[[1]] > 0 & my.Trait.Cor[[2]] > 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = positive, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = positive, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

    
# MAKE PLOT (META-ANALYSIS)

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(consensusCor, 2), "\n(",
                        signif(consensusPvalue, 1), ")", sep = "")

    dim(textMatrix) = dim(consensusCor)

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Transposition_Regulator_SNVs_META.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = consensusCor,
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with L1 Regulators SNVs (EUR & YRI)"))
    
    dev.off()


    
    
                 
# CORRELATION TO POLYMORPHIC SVs --------------------------------------------
     
    

# Define the variants of interest
my.variants <- c('INV_delly_INV00066128', 'INV_delly_INV00003623', 'ALU_umary_ALU_7919', 'ALU_umary_ALU_3176', 'ALU_umary_ALU_8971', 'L1_umary_LINE1_1066')

# Define plot labels for variants
my.variant.labels <- c('INV_delly_INV00066128', 'INV_delly_INV00003623', 'ALU_umary_ALU_7919', 'ALU_umary_ALU_3176', 'ALU_umary_ALU_8971', 'L1_umary_LINE1_1066') 

# Extract genotypes and define them as the phenotype
my.phenotype.EUR <- as.data.frame(binary_genotypes.EUR[rownames(my.consMEs[[1]][["data"]]), my.variants])
my.phenotype.YRI <- as.data.frame(binary_genotypes.YRI[rownames(my.consMEs[[2]][["data"]]), my.variants])

# Put the phenotype data into a multiset 
my.Traits = vector(mode = "list", length = nSets)

    my.Traits[[1]] = list(data = as.data.frame(my.phenotype.EUR))
    names(my.Traits[[1]]$data) = colnames(my.phenotype.EUR)
    rownames(my.Traits[[1]]$data) = rownames(my.phenotype.EUR)
    
    my.Traits[[2]] = list(data = as.data.frame(my.phenotype.YRI))
    names(my.Traits[[2]]$data) = colnames(my.phenotype.YRI)
    rownames(my.Traits[[2]]$data) = rownames(my.phenotype.YRI)

# Set up variables to contain the module-trait correlations/pvalues
my.Trait.Cor = list()
my.Trait.Pvalue = list()

# Calculate the correlations/pvalues
for (set in 1:nSets){
  
    my.Trait.Cor[[set]] = cor(my.consMEs[[set]]$data, my.Traits[[set]]$data, use = "p")
    
    my.Trait.Pvalue[[set]] = corPvalueFisher(my.Trait.Cor[[set]], nrow(my.consMEs[[set]][["data"]]))
    
}


# MAKE PLOT (SET1)

    # Define the set number
    set = 1

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Polymorphic_SVs_EUR.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Polymorphic SVs in", setLabels[set]))
    
    dev.off()

    
# MAKE PLOT (SET2)

    # Define the set number
    set = 2

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(my.Trait.Cor[[set]], 2), "\n(",
                        signif(my.Trait.Pvalue[[set]], 1), ")", sep = "")

    dim(textMatrix) = dim(my.Trait.Cor[[set]])

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Polymorphic_SVs_YRI.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = my.Trait.Cor[[set]],
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Polymorphic SVs in", setLabels[set]))
    
    dev.off()





# COMBINE RESULTS
    
# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))
consensusPvalue = matrix(NA, nrow(my.Trait.Cor[[1]]), ncol(my.Trait.Cor[[1]]))

# Find consensus negative correlations
negative = my.Trait.Cor[[1]] < 0 & my.Trait.Cor[[2]] < 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = negative, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = negative, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

# Find consensus positive correlations
positive = my.Trait.Cor[[1]] > 0 & my.Trait.Cor[[2]] > 0

consensusCor <- trait_cor_meta(consistent.logical.matrix = positive, 
                               input.consensus = consensusCor, 
                               input.multiset = my.Trait.Cor, 
                               multiset_type = 'correlation')

consensusPvalue <- trait_cor_meta(consistent.logical.matrix = positive, 
                                  input.consensus = consensusPvalue, 
                                  input.multiset = my.Trait.Pvalue, 
                                  multiset_type = 'pvalue')

    
# MAKE PLOT (META-ANALYSIS)

    # Define text matrix for correlations/pvalues
    textMatrix =  paste(signif(consensusCor, 2), "\n(",
                        signif(consensusPvalue, 1), ")", sep = "")

    dim(textMatrix) = dim(consensusCor)

    # plot results
    pdf(file = paste(dir.trait.cor, "Consensus_Net_vs_Polymorphic_SVs_META.pdf", sep = ''), wi = 10, he = 14);
    par(mar = c(6, 8.8, 3, 2.2))
    
    labeledHeatmap(Matrix = consensusCor,
                   xLabels = names(my.Traits[[set]]$data),
                   yLabels = MEColorNames,
                   ySymbols = MEColorNames,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(300),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 1.25,
                   zlim = c(-1,1),
                   main = paste("Correlations with Polymorphic SVs (EUR & YRI)"))
    
    dev.off()


    
    
                 

# Session Info ------------------------------------------------------------------------------------------------------------------------



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_WGCNA_Trait_Correlations.txt", sep =""))
sessionInfo()
sink()      
    

    
# Clean the environment
rm(list=ls())


