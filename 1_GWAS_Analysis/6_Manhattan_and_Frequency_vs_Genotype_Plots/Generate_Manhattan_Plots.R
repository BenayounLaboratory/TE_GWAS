# Set strings as factors
options(stringsAsFactors = F)

# Load libraries  
library(data.table) # for fread and fsave
library(CMplot) # Manhattan plots

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/6_Manhattan_and_Frequency_vs_Genotype_Plots/Generate_Manhattan_Plots_Functions.R")     

# Define output directory
output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/6_Manhattan_and_Frequency_vs_Genotype_Plots/Manhattan_Plots/'





# DEFINE INPUT/GENERAL PARAMETERS

# Load table with empirical pval/FDR thresholds
empirical.thresholds <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/EMPIRICAL_PVALUE_AND_FDR_THRESHOLDS.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

# Load ALL BACKGROUND SNV and SV annotations
background.SNV <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/All_Annotated_GWAS_Background_SNVs_only.csv", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = ',')
background.SV <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/All_Annotated_GWAS_Background_SVs_only.csv", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = ',')
background.all <- rbind(background.SNV, background.SV)

# Define the background thats in the ENCODE blacklist
background.blacklist <- background.all[!is.na(background.all$Blacklist_Label),]

# Define FDR thresholds
BH.FDR.threshold <- 0.05

    



# Plot: L1+Alu combined singleton GWAS (META)

# Read external files. Only keep stats.
Meta.stats <- fread('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/plink.meta', header = TRUE, sep = ' ', data.table = FALSE)

    # Add BH FDR to stats
    Meta.stats$FDR <- p.adjust(Meta.stats$`P(R)`, method = 'BH')
    
    # Only plot results with pvalue smaller than 1e-2 (to reduce plotting time and filesize)
    Meta.stats <- Meta.stats[which(Meta.stats$`P(R)` <= 1e-2), ]
    
    # Update the row numbers
    rownames(Meta.stats) <- 1:nrow(Meta.stats)

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(Meta.stats[Meta.stats$FDR < BH.FDR.threshold, 'P(R)']) 
BH.pval_threshold # 5.995e-06

# Define pvalue at average empirical FDR threshold
empirical.pval_threshold <- empirical.thresholds[c('Combined_meta'), c('P_threshold')]
empirical.pval_threshold # 1.396e-05

# Transform to CMplot format
CMplot_GWAS <- Meta_to_CMplot(meta.results = Meta.stats)

# Load blacklist variants (to label on the Manhattan plot)
blacklist.snps <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/Annotated_Significant_Variants_Blacklist_SNVs.csv", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = ',')

    # Find CMPlot rows for significant SNVs
    sig.snp.rows <- which(CMplot_GWAS$SNP %in% blacklist.snps$SNP)
    
    # Map CMPlot rows to CMPlot indices
    sig.snp.indices <- CMplot_GWAS[sig.snp.rows, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Manhattan_Plot_Combined_Singleton_Meta-Analysis_FDR5", '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_GWAS[, -c(5)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           axis.lwd = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 5),
           threshold.col = c('red', 'red'),
           ylim = c(2, 30),
           highlight = sig.snp.indices,
           highlight.col = c('blue'),
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           #highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           lab.cex = 1,
           axis.cex = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()


# Loop over permutation files and plot

# Define output directory
perm.output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/6_Manhattan_and_Frequency_vs_Genotype_Plots/Manhattan_Plots/Combined_Singleton_GWAS_Permutations/'

# Define permutation data files
permutations.dir <- c("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations")
permutations.list <- list.files(permutations.dir, '\\.meta$', recursive = TRUE, full.names = TRUE)

# Define vector with permutation number (note that files are loaded alphabetically and not from smallest to largest number)
permutation_labels <- sort(as.character(1:length(permutations.list)))

# Define loop counter
counter <- 0

for (permutation_file in permutations.list) {
  
    # Update counter 
    counter <- counter + 1
    
    # Extract current permutation label
    ith_permutation_label <- permutation_labels[counter]
  
    # Read external files. Only keep stats.
    #Meta.stats <- read.csv(permutation_file, header = TRUE, sep = ' ')
    Meta.stats <- fread(permutation_file, header = TRUE, sep = ' ', data.table = FALSE)

    # Only plot results with pvalue smaller than 1e-1 (to reduce plotting time and filesize)
    Meta.stats <- Meta.stats[which(Meta.stats$`P(R)` <= 1e-2), ]
    
    # Transform to CMplot format
    CMplot_GWAS <- Meta_to_CMplot(meta.results = Meta.stats)
    
    # Make plots and save to file
    pdf(file = paste(perm.output.dir, "Manhattan_Plot_Combined_Singleton_Meta-Analysis_FDR5_PERMUTATION_", ith_permutation_label, '.pdf', sep=""), width = 10, height = 5)
    par(mar = c(5.1,6.1,4.1,2.1))
        CMplot(CMplot_GWAS[, -c(5)],
               band = 0,
               plot.type = "m",
               type = "p",
               col = c("black","grey60"),
               LOG10 = TRUE,
               amplify = FALSE,
               axis.lwd = 0.75,
               threshold.lwd = 0.75,
               threshold = c(BH.pval_threshold, empirical.pval_threshold),
               threshold.lty = c(1, 5),
               threshold.col = c('red', 'red'),
               ylim = c(2, 30),
               #highlight = highlight_indices[1:20],
               #highlight.col = NULL,
               #highlight.cex = 0.25,
               #highlight.text = highlight_genes[1:20],
               #highlight.text.cex = 0.75,
               cex = 0.25,
               lab.cex = 1,
               axis.cex = 1,
               chr.labels.angle = 45,
               #width = 20,
               #height = 10,
               #file = "tiff",
               verbose = TRUE,
               #dpi = 300,
               file.output = F)
    dev.off()
}

        
      
    






# Plot: L1+Alu combined singleton GWAS (AFR)

# Read external files. Only keep stats.
gwas.stats <- fread('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AFR_GWAS.assoc.logistic', header = TRUE, sep = ' ', data.table = FALSE)

    # Add BH FDR to stats
    gwas.stats$FDR <- p.adjust(gwas.stats$P, method = 'BH')
    
    # Only plot results with pvalue smaller than 1e-2 (to reduce plotting time and filesize)
    gwas.stats <- gwas.stats[which(gwas.stats$P <= 1e-2), ]
    
    # Update the row numbers
    rownames(gwas.stats) <- 1:nrow(gwas.stats)

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(gwas.stats[gwas.stats$FDR < BH.FDR.threshold, 'P']) 
BH.pval_threshold # 1.177e-06

# Define pvalue at average empirical FDR threshold
empirical.pval_threshold <- empirical.thresholds[c('Combined_AFR'), c('P_threshold')]
empirical.pval_threshold # 4.606e-06

# Define the stricter threshold
stricter.threshold <- min(BH.pval_threshold, empirical.pval_threshold)

# Transform to CMplot format
CMplot_GWAS <- Individual_to_CMplot(gwas.results = gwas.stats)

# color-code blacklisted variants on the Manhattan plot

    # Find significant variants that are also blacklisted
    sig.and.blacklisted <- intersect(CMplot_GWAS[CMplot_GWAS$P <= stricter.threshold, 'SNP'], background.blacklist$RsID)
    
    # Find CMPlot rows for significant blacklisted variants
    sig.snp.rows <- which(CMplot_GWAS$SNP %in% sig.and.blacklisted)
    
    # Map CMPlot rows to CMPlot indices
    sig.snp.indices <- CMplot_GWAS[sig.snp.rows, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Manhattan_Plot_Combined_Singleton_GWAS_AFR-only_FDR5", '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_GWAS[, -c(5)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           axis.lwd = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 5),
           threshold.col = c('red', 'red'),
           ylim = c(2, 20),
           highlight = sig.snp.indices,
           highlight.col = c('blue'),
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           #highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           lab.cex = 1,
           axis.cex = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()

        
      
    






# Plot: L1+Alu combined singleton GWAS (AMR)

# Read external files. Only keep stats.
gwas.stats <- fread('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AMR_GWAS.assoc.logistic', header = TRUE, sep = ' ', data.table = FALSE)

    # Add BH FDR to stats
    gwas.stats$FDR <- p.adjust(gwas.stats$P, method = 'BH')
    
    # Only plot results with pvalue smaller than 1e-2 (to reduce plotting time and filesize)
    gwas.stats <- gwas.stats[which(gwas.stats$P <= 1e-2), ]
    
    # Update the row numbers
    rownames(gwas.stats) <- 1:nrow(gwas.stats)

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(gwas.stats[gwas.stats$FDR < BH.FDR.threshold, 'P']) 
BH.pval_threshold # 3.526e-07

# Define pvalue at average empirical FDR threshold
empirical.pval_threshold <- empirical.thresholds[c('Combined_AMR'), c('P_threshold')]
empirical.pval_threshold # 1.069e-06

# Define the stricter threshold
stricter.threshold <- min(BH.pval_threshold, empirical.pval_threshold)

# Transform to CMplot format
CMplot_GWAS <- Individual_to_CMplot(gwas.results = gwas.stats)

# color-code blacklisted variants on the Manhattan plot

    # Find significant variants that are also blacklisted
    sig.and.blacklisted <- intersect(CMplot_GWAS[CMplot_GWAS$P <= stricter.threshold, 'SNP'], background.blacklist$RsID)
    
    # Find CMPlot rows for significant blacklisted variants
    sig.snp.rows <- which(CMplot_GWAS$SNP %in% sig.and.blacklisted)
    
    # Map CMPlot rows to CMPlot indices
    sig.snp.indices <- CMplot_GWAS[sig.snp.rows, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Manhattan_Plot_Combined_Singleton_GWAS_AMR-only_FDR5", '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_GWAS[, -c(5)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           axis.lwd = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 5),
           threshold.col = c('red', 'red'),
           ylim = c(2, 20),
           highlight = sig.snp.indices,
           highlight.col = c('blue'),
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           #highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           lab.cex = 1,
           axis.cex = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()

        
      
    






# Plot: L1+Alu combined singleton GWAS (EAS)

# Read external files. Only keep stats.
gwas.stats <- fread('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EAS_GWAS.assoc.logistic', header = TRUE, sep = ' ', data.table = FALSE)

    # Add BH FDR to stats
    gwas.stats$FDR <- p.adjust(gwas.stats$P, method = 'BH')
    
    # Only plot results with pvalue smaller than 1e-2 (to reduce plotting time and filesize)
    gwas.stats <- gwas.stats[which(gwas.stats$P <= 1e-2), ]
    
    # Update the row numbers
    rownames(gwas.stats) <- 1:nrow(gwas.stats)

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(gwas.stats[gwas.stats$FDR < BH.FDR.threshold, 'P']) 
BH.pval_threshold # NA

# Define pvalue at average empirical FDR threshold
empirical.pval_threshold <- empirical.thresholds[c('Combined_EAS'), c('P_threshold')]
empirical.pval_threshold # NA

# Define the stricter threshold
stricter.threshold <- min(BH.pval_threshold, empirical.pval_threshold)

# Transform to CMplot format
CMplot_GWAS <- Individual_to_CMplot(gwas.results = gwas.stats)

# color-code blacklisted variants on the Manhattan plot

    # Find significant variants that are also blacklisted
    sig.and.blacklisted <- intersect(CMplot_GWAS[CMplot_GWAS$P <= stricter.threshold, 'SNP'], background.blacklist$RsID)
    
    # Find CMPlot rows for significant blacklisted variants
    sig.snp.rows <- which(CMplot_GWAS$SNP %in% sig.and.blacklisted)
    
    # Map CMPlot rows to CMPlot indices
    sig.snp.indices <- CMplot_GWAS[sig.snp.rows, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Manhattan_Plot_Combined_Singleton_GWAS_EAS-only_FDR5", '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_GWAS[, -c(5)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           axis.lwd = 0.75,
           threshold.lwd = 0.75,
           #threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 5),
           threshold.col = c('red', 'red'),
           ylim = c(2, 20),
           highlight = sig.snp.indices,
           highlight.col = c('blue'),
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           #highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           lab.cex = 1,
           axis.cex = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()

        
      
    






# Plot: L1+Alu combined singleton GWAS (EUR)

# Read external files. Only keep stats.
gwas.stats <- fread('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EUR_GWAS.assoc.logistic', header = TRUE, sep = ' ', data.table = FALSE)

    # Add BH FDR to stats
    gwas.stats$FDR <- p.adjust(gwas.stats$P, method = 'BH')
    
    # Only plot results with pvalue smaller than 1e-2 (to reduce plotting time and filesize)
    gwas.stats <- gwas.stats[which(gwas.stats$P <= 1e-2), ]
    
    # Update the row numbers
    rownames(gwas.stats) <- 1:nrow(gwas.stats)

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(gwas.stats[gwas.stats$FDR < BH.FDR.threshold, 'P']) 
BH.pval_threshold # NA

# Define pvalue at average empirical FDR threshold
empirical.pval_threshold <- empirical.thresholds[c('Combined_EUR'), c('P_threshold')]
empirical.pval_threshold # NA

# Define the stricter threshold
stricter.threshold <- min(BH.pval_threshold, empirical.pval_threshold)

# Transform to CMplot format
CMplot_GWAS <- Individual_to_CMplot(gwas.results = gwas.stats)

# color-code blacklisted variants on the Manhattan plot

    # Find significant variants that are also blacklisted
    sig.and.blacklisted <- intersect(CMplot_GWAS[CMplot_GWAS$P <= stricter.threshold, 'SNP'], background.blacklist$RsID)
    
    # Find CMPlot rows for significant blacklisted variants
    sig.snp.rows <- which(CMplot_GWAS$SNP %in% sig.and.blacklisted)
    
    # Map CMPlot rows to CMPlot indices
    sig.snp.indices <- CMplot_GWAS[sig.snp.rows, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Manhattan_Plot_Combined_Singleton_GWAS_EUR-only_FDR5", '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_GWAS[, -c(5)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           axis.lwd = 0.75,
           threshold.lwd = 0.75,
           #threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 5),
           threshold.col = c('red', 'red'),
           ylim = c(2, 20),
           highlight = sig.snp.indices,
           highlight.col = c('blue'),
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           #highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           lab.cex = 1,
           axis.cex = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()

        
      
    






# Plot: L1+Alu combined singleton GWAS (SAS)

# Read external files. Only keep stats.
gwas.stats <- fread('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/SAS_GWAS.assoc.logistic', header = TRUE, sep = ' ', data.table = FALSE)

    # Add BH FDR to stats
    gwas.stats$FDR <- p.adjust(gwas.stats$P, method = 'BH')
    
    # Only plot results with pvalue smaller than 1e-2 (to reduce plotting time and filesize)
    gwas.stats <- gwas.stats[which(gwas.stats$P <= 1e-2), ]
    
    # Update the row numbers
    rownames(gwas.stats) <- 1:nrow(gwas.stats)

# Extract pvalue corresponding to the desired BH FDR
BH.pval_threshold <- max(gwas.stats[gwas.stats$FDR < BH.FDR.threshold, 'P']) 
BH.pval_threshold # 9.056e-07

# Define pvalue at average empirical FDR threshold
empirical.pval_threshold <- empirical.thresholds[c('Combined_SAS'), c('P_threshold')]
empirical.pval_threshold # 8.463e-07

# Define the stricter threshold
stricter.threshold <- min(BH.pval_threshold, empirical.pval_threshold)

# Transform to CMplot format
CMplot_GWAS <- Individual_to_CMplot(gwas.results = gwas.stats)

# color-code blacklisted variants on the Manhattan plot

    # Find significant variants that are also blacklisted
    sig.and.blacklisted <- intersect(CMplot_GWAS[CMplot_GWAS$P <= stricter.threshold, 'SNP'], background.blacklist$RsID)
    
    # Find CMPlot rows for significant blacklisted variants
    sig.snp.rows <- which(CMplot_GWAS$SNP %in% sig.and.blacklisted)
    
    # Map CMPlot rows to CMPlot indices
    sig.snp.indices <- CMplot_GWAS[sig.snp.rows, 'index']
    
# Make plot and save to file
pdf(file = paste(output.dir, "Manhattan_Plot_Combined_Singleton_GWAS_SAS-only_FDR5", '.pdf', sep=""), width = 10, height = 5)
par(mar = c(5.1,6.1,4.1,2.1)) 

    CMplot(CMplot_GWAS[, -c(5)],
           band = 0,
           plot.type = "m", 
           type = "p",
           col = c("black","grey60"),
           LOG10 = TRUE,
           amplify = FALSE,
           axis.lwd = 0.75,
           threshold.lwd = 0.75,
           threshold = c(BH.pval_threshold, empirical.pval_threshold),
           threshold.lty = c(1, 5),
           threshold.col = c('red', 'red'),
           ylim = c(2, 20),
           highlight = sig.snp.indices,
           highlight.col = c('blue'),
           highlight.cex = 0.25,
           #highlight.text = aggregate.annotation$symbol,
           #highlight.text.cex = 0.60,
           #highlight.text.xadj = c(1, -1, -1, -1, 1, -1, -1),
           #highlight.text.yadj = c(1, 0, 0, 1, 0, 1, 0),
           cex = 0.25, 
           lab.cex = 1,
           axis.cex = 1,
           chr.labels.angle = 45,
           #width = 20,
           #height = 10,
           #file = "tiff",
           verbose = TRUE,
           #dpi = 300,
           file.output = F)
dev.off()





# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/6_Manhattan_and_Frequency_vs_Genotype_Plots/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Manhattan_Plots.txt", sep =""))
sessionInfo()
sink()      



# Clean the environment
rm(list=ls())       

