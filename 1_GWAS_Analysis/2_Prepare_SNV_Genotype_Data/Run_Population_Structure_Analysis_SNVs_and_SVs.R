# Set strings as factors
options(stringsAsFactors = F)





# COMBINED SINGLE NUCLEOTIDE VARIANT (SNV) AND STUCTURAL VARIANT (SV) POPULATION STRUCTURE ANALYSIS USING PLINK 

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Population_Structure_Analysis/'

# Load SV metadata 
SV_sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/integrated_call_samples_v3.20130502.ALL.panel.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(SV_sample_info) <- SV_sample_info$sample

    # Remove samples that do not have both SNV and SV data
    SV_sample_info <- SV_sample_info[-which(SV_sample_info$sample %in% c('NA18498')), ]

# Extract metadata by superpopulation
AFR.samples <- SV_sample_info[which(SV_sample_info$super_pop == 'AFR'), ]
AMR.samples <- SV_sample_info[which(SV_sample_info$super_pop == 'AMR'), ]
EAS.samples <- SV_sample_info[which(SV_sample_info$super_pop == 'EAS'), ]
EUR.samples <- SV_sample_info[which(SV_sample_info$super_pop == 'EUR'), ]
SAS.samples <- SV_sample_info[which(SV_sample_info$super_pop == 'SAS'), ]

# Load PCA eigenvectors.
pca.AFR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_PCA/Final_SNV_SV_AFR_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")
pca.AMR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_PCA/Final_SNV_SV_AMR_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")
pca.EAS <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_PCA/Final_SNV_SV_EAS_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")
pca.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_PCA/Final_SNV_SV_EUR_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")
pca.SAS <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_PCA/Final_SNV_SV_SAS_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")

    # Match the order of the sample names
    pca.AFR <- pca.AFR[rownames(AFR.samples), ]
    pca.AMR <- pca.AMR[rownames(AMR.samples), ]
    pca.EAS <- pca.EAS[rownames(EAS.samples), ]
    pca.EUR <- pca.EUR[rownames(EUR.samples), ]
    pca.SAS <- pca.SAS[rownames(SAS.samples), ]
    
    # Clean PCA eigenvector df: remove sample names column 
    pca.AFR <- pca.AFR[, -c(1)]
    pca.AMR <- pca.AMR[, -c(1)]
    pca.EAS <- pca.EAS[, -c(1)]
    pca.EUR <- pca.EUR[, -c(1)]
    pca.SAS <- pca.SAS[, -c(1)]
    
    # Clean PCA eigenvector df: assign PC labels
    colnames(pca.AFR) <- paste('PC', 1:ncol(pca.EUR), sep = '')
    colnames(pca.AMR) <- paste('PC', 1:ncol(pca.AMR), sep = '')
    colnames(pca.EAS) <- paste('PC', 1:ncol(pca.EUR), sep = '')
    colnames(pca.EUR) <- paste('PC', 1:ncol(pca.EUR), sep = '')
    colnames(pca.SAS) <- paste('PC', 1:ncol(pca.SAS), sep = '')
    
# Load PCA eigenvalues.
eigenval.AFR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_PCA/Final_SNV_SV_AFR_PCA.eigenval", header = FALSE)
eigenval.AMR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_PCA/Final_SNV_SV_AMR_PCA.eigenval", header = FALSE)
eigenval.EAS <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_PCA/Final_SNV_SV_EAS_PCA.eigenval", header = FALSE)
eigenval.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_PCA/Final_SNV_SV_EUR_PCA.eigenval", header = FALSE)
eigenval.SAS <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_PCA/Final_SNV_SV_SAS_PCA.eigenval", header = FALSE)

# Combine PCs and save
pca.ALL <- rbind(pca.AFR, pca.AMR, pca.EAS, pca.EUR, pca.SAS)
write.table(pca.ALL, file = paste(dir.output, "SNV_AND_SV_Population_Structure_PCA", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)

# Calculate percent variance explained (pve)
pve.AFR <- 100 * (eigenval.AFR/sum(eigenval.AFR)) 
pve.AMR <- 100 * (eigenval.AMR/sum(eigenval.AMR)) 
pve.EAS <- 100 * (eigenval.EAS/sum(eigenval.EAS)) 
pve.EUR <- 100 * (eigenval.EUR/sum(eigenval.EUR)) 
pve.SAS <- 100 * (eigenval.SAS/sum(eigenval.SAS)) 

# Add a 'color' column to the sample tables
AFR.samples$Point_color <- '#648FFF'
AMR.samples$Point_color <- '#785EF0'
EAS.samples$Point_color <- '#DC267F'
EUR.samples$Point_color <- '#FE6100'
SAS.samples$Point_color <- '#FFB000'

    # AFR: Update the point colors
    AFR.samples[which(AFR.samples$pop == 'ACB'), 'Point_color'] <- '#0A0E19'
    AFR.samples[which(AFR.samples$pop == 'GWD'), 'Point_color'] <- '#1E2B4D'
    AFR.samples[which(AFR.samples$pop == 'ESN'), 'Point_color'] <- '#324880'
    AFR.samples[which(AFR.samples$pop == 'MSL'), 'Point_color'] <- '#4664B3'
    AFR.samples[which(AFR.samples$pop == 'YRI'), 'Point_color'] <- '#5A81E6'
    AFR.samples[which(AFR.samples$pop == 'LWK'), 'Point_color'] <- '#7BA0FF'
    AFR.samples[which(AFR.samples$pop == 'ASW'), 'Point_color'] <- '#9AB6FF'
    
    # AMR: Update the point colors
    AMR.samples[which(AMR.samples$pop == 'PUR'), 'Point_color'] <- '#0C0918'
    AMR.samples[which(AMR.samples$pop == 'CLM'), 'Point_color'] <- '#423484'
    AMR.samples[which(AMR.samples$pop == 'PEL'), 'Point_color'] <- '#7F66F1'
    AMR.samples[which(AMR.samples$pop == 'MXL'), 'Point_color'] <- '#BCAFF8'

    # EAS: Update the point colors
    EAS.samples[which(EAS.samples$pop == 'CHS'), 'Point_color'] <- '#16040D'
    EAS.samples[which(EAS.samples$pop == 'CDX'), 'Point_color'] <- '#6E1340'
    EAS.samples[which(EAS.samples$pop == 'KHV'), 'Point_color'] <- '#C62272'
    EAS.samples[which(EAS.samples$pop == 'CHB'), 'Point_color'] <- '#E872AC'
    EAS.samples[which(EAS.samples$pop == 'JPT'), 'Point_color'] <- '#F6C9DF'
    
    # EUR: Update the point colors
    EUR.samples[which(EUR.samples$pop == 'GBR'), 'Point_color'] <- '#190A00'
    EUR.samples[which(EUR.samples$pop == 'FIN'), 'Point_color'] <- '#7F3100'
    EUR.samples[which(EUR.samples$pop == 'IBS'), 'Point_color'] <- '#E55700'
    EUR.samples[which(EUR.samples$pop == 'CEU'), 'Point_color'] <- '#FE9859'
    EUR.samples[which(EUR.samples$pop == 'TSI'), 'Point_color'] <- '#FFD8BF'
    
    # SAS: Update the point colors
    SAS.samples[which(SAS.samples$pop == 'PJL'), 'Point_color'] <- '#191200'
    SAS.samples[which(SAS.samples$pop == 'BEB'), 'Point_color'] <- '#805800'
    SAS.samples[which(SAS.samples$pop == 'STU'), 'Point_color'] <- '#E69E00'
    SAS.samples[which(SAS.samples$pop == 'ITU'), 'Point_color'] <- '#FFCC59'
    SAS.samples[which(SAS.samples$pop == 'GIH'), 'Point_color'] <- '#FFE3A6'

# Add a point 'shape' column to the sample table
AFR.samples$Point_pch <- as.numeric('16')
AMR.samples$Point_pch <- as.numeric('16')
EAS.samples$Point_pch <- as.numeric('16')
EUR.samples$Point_pch <- as.numeric('16')
SAS.samples$Point_pch <- as.numeric('16')
    
    # AFR: Update the point shapes
    AFR.samples[which(AFR.samples$pop == 'ACB'), 'Point_pch'] <- 15
    AFR.samples[which(AFR.samples$pop == 'GWD'), 'Point_pch'] <- 16
    AFR.samples[which(AFR.samples$pop == 'ESN'), 'Point_pch'] <- 17
    AFR.samples[which(AFR.samples$pop == 'MSL'), 'Point_pch'] <- 18
    AFR.samples[which(AFR.samples$pop == 'YRI'), 'Point_pch'] <- 3
    AFR.samples[which(AFR.samples$pop == 'LWK'), 'Point_pch'] <- 8
    AFR.samples[which(AFR.samples$pop == 'ASW'), 'Point_pch'] <- 9
    
    # AMR: Update the point shapes
    AMR.samples[which(AMR.samples$pop == 'PUR'), 'Point_pch'] <- 15
    AMR.samples[which(AMR.samples$pop == 'CLM'), 'Point_pch'] <- 16
    AMR.samples[which(AMR.samples$pop == 'PEL'), 'Point_pch'] <- 17
    AMR.samples[which(AMR.samples$pop == 'MXL'), 'Point_pch'] <- 18

    # EAS: Update the point shapes
    EAS.samples[which(EAS.samples$pop == 'CHS'), 'Point_pch'] <- 15
    EAS.samples[which(EAS.samples$pop == 'CDX'), 'Point_pch'] <- 16
    EAS.samples[which(EAS.samples$pop == 'KHV'), 'Point_pch'] <- 17
    EAS.samples[which(EAS.samples$pop == 'CHB'), 'Point_pch'] <- 18
    EAS.samples[which(EAS.samples$pop == 'JPT'), 'Point_pch'] <- 3
    
    # EUR: Update the point shapes
    EUR.samples[which(EUR.samples$pop == 'GBR'), 'Point_pch'] <- 15
    EUR.samples[which(EUR.samples$pop == 'FIN'), 'Point_pch'] <- 16
    EUR.samples[which(EUR.samples$pop == 'IBS'), 'Point_pch'] <- 17
    EUR.samples[which(EUR.samples$pop == 'CEU'), 'Point_pch'] <- 18
    EUR.samples[which(EUR.samples$pop == 'TSI'), 'Point_pch'] <- 3
    
    # SAS: Update the point shapes
    SAS.samples[which(SAS.samples$pop == 'PJL'), 'Point_pch'] <- 15
    SAS.samples[which(SAS.samples$pop == 'BEB'), 'Point_pch'] <- 16
    SAS.samples[which(SAS.samples$pop == 'STU'), 'Point_pch'] <- 17
    SAS.samples[which(SAS.samples$pop == 'ITU'), 'Point_pch'] <- 18
    SAS.samples[which(SAS.samples$pop == 'GIH'), 'Point_pch'] <- 3

# Define plot and legend parameters
AFR.labels <- c('ACB', 'GWD', 'ESN', 'MSL', 'YRI', 'LWK', 'ASW')
AFR.colors <- c('#0A0E19', '#1E2B4D', '#324880',  '#4664B3', '#5A81E6', '#7BA0FF', '#9AB6FF')
AFR.pch <- c(15, 16, 17, 18, 3, 8, 9)

AMR.labels <- c('PUR', 'CLM', 'PEL', 'MXL')
AMR.colors <- c('#0C0918', '#423484', '#7F66F1', '#BCAFF8')
AMR.pch <- c(15, 16, 17, 18)

EAS.labels <- c('CHS', 'CDX', 'KHV', 'CHB', 'JPT')
EAS.colors <- c('#16040D', '#6E1340', '#C62272', '#E872AC', '#F6C9DF')
EAS.pch <- c(15, 16, 17, 18, 3)

EUR.labels <- c('GBR', 'FIN', 'IBS', 'CEU', 'TSI')
EUR.colors <- c('#190A00', '#7F3100', '#E55700', '#FE9859', '#FFD8BF')
EUR.pch <- c(15, 16, 17, 18, 3)

SAS.labels <- c('PJL', 'BEB', 'STU', 'ITU', 'GIH')
SAS.colors <- c('#191200', '#805800', '#E69E00', '#FFCC59', '#FFE3A6')
SAS.pch <- c(15, 16, 17, 18, 3)

pca.label.cex <- 1
pca.axis.cex <- 1


# Plot pve

    # AFR
    pdf(paste(dir.output, "Barplot_SNV_AND_SV_PCA_percent_variance_AFR", ".pdf", sep=""))
    par(mfrow=c(2,1))
        barplot(t(pve.AFR), names.arg = rownames(pve.AFR), xlab = "Principal component", ylim = c(0, 20), ylab = 'Percent variance explained', main = 'AFR SNV+SV PCA', las = 2)
    dev.off()
    
    # AMR
    pdf(paste(dir.output, "Barplot_SNV_AND_SV_PCA_percent_variance_AMR", ".pdf", sep=""))
    par(mfrow=c(2,1))
        barplot(t(pve.AMR), names.arg = rownames(pve.AMR), xlab = "Principal component", ylim = c(0, 35), ylab = 'Percent variance explained', main = 'AMR SNV+SV PCA', las = 2)
    dev.off()
    
    # EAS
    pdf(paste(dir.output, "Barplot_SNV_AND_SV_PCA_percent_variance_EAS", ".pdf", sep=""))
    par(mfrow=c(2,1))
        barplot(t(pve.EAS), names.arg = rownames(pve.EAS), xlab = "Principal component", ylim = c(0, 20), ylab = 'Percent variance explained', main = 'EAS SNV+SV PCA', las = 2)
    dev.off()
    
    # EUR
    pdf(paste(dir.output, "Barplot_SNV_AND_SV_PCA_percent_variance_EUR", ".pdf", sep=""))
    par(mfrow=c(2,1))
        barplot(t(pve.EUR), names.arg = rownames(pve.EUR), xlab = "Principal component", ylim = c(0, 15), ylab = 'Percent variance explained', main = 'EUR SNV+SV PCA', las = 2)
    dev.off()
    
    # SAS
    pdf(paste(dir.output, "Barplot_SNV_AND_SV_PCA_percent_variance_SAS", ".pdf", sep=""))
    par(mfrow=c(2,1))
        barplot(t(pve.SAS), names.arg = rownames(pve.SAS), xlab = "Principal component", ylim = c(0, 10), ylab = 'Percent variance explained', main = 'SAS SNV+SV PCA', las = 2)
    dev.off()


# Plot PCs vs Ancestry
    
    # AFR
    pdf(paste(dir.output, "Plot_SNV_AND_SV_PCs_vs_Ancestry_AFR",".pdf", sep=""), width = 9, height = 9)
    par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses
    
        # PC1 vs PC2
        plot(pca.AFR$PC1,
             pca.AFR$PC2,
             col = AFR.samples$Point_color,
             pch = AFR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC2', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
            # add legend
            legend('bottomleft', legend = AFR.labels, col = AFR.colors, pch = AFR.pch, cex = 1, lty = 0, bty = 'n')
        
        # PC1 vs PC3
        plot(pca.AFR$PC1,
             pca.AFR$PC3,
             col = AFR.samples$Point_color,
             pch = AFR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC3', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC4
        plot(pca.AFR$PC1,
             pca.AFR$PC4,
             col = AFR.samples$Point_color,
             pch = AFR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC4', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC5
        plot(pca.AFR$PC1,
             pca.AFR$PC5,
             col = AFR.samples$Point_color,
             pch = AFR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC5', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
    dev.off()
    
    # AMR
    pdf(paste(dir.output, "Plot_SNV_AND_SV_PCs_vs_Ancestry_AMR",".pdf", sep=""), width = 9, height = 9)
    par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses
    
        # PC1 vs PC2
        plot(pca.AMR$PC1,
             pca.AMR$PC2,
             col = AMR.samples$Point_color,
             pch = AMR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC2', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
            # add legend
            legend('topleft', legend = AMR.labels, col = AMR.colors, pch = AMR.pch, cex = 1, lty = 0, bty = 'n')
        
        # PC1 vs PC3
        plot(pca.AMR$PC1,
             pca.AMR$PC3,
             col = AMR.samples$Point_color,
             pch = AMR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC3', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC4
        plot(pca.AMR$PC1,
             pca.AMR$PC4,
             col = AMR.samples$Point_color,
             pch = AMR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC4', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC5
        plot(pca.AMR$PC1,
             pca.AMR$PC5,
             col = AMR.samples$Point_color,
             pch = AMR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC5', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
    dev.off()
    
    # EAS
    pdf(paste(dir.output, "Plot_SNV_AND_SV_PCs_vs_Ancestry_EAS",".pdf", sep=""), width = 9, height = 9)
    par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses
    
        # PC1 vs PC2
        plot(pca.EAS$PC1,
             pca.EAS$PC2,
             col = EAS.samples$Point_color,
             pch = EAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC2', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
            # add legend
            legend('topleft', legend = EAS.labels, col = EAS.colors, pch = EAS.pch, cex = 1, lty = 0, bty = 'n')
        
        # PC1 vs PC3
        plot(pca.EAS$PC1,
             pca.EAS$PC3,
             col = EAS.samples$Point_color,
             pch = EAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC3', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC4
        plot(pca.EAS$PC1,
             pca.EAS$PC4,
             col = EAS.samples$Point_color,
             pch = EAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC4', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC5
        plot(pca.EAS$PC1,
             pca.EAS$PC5,
             col = EAS.samples$Point_color,
             pch = EAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC5', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
    dev.off()
    
    # EUR
    pdf(paste(dir.output, "Plot_SNV_AND_SV_PCs_vs_Ancestry_EUR",".pdf", sep=""), width = 9, height = 9)
    par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses
    
        # PC1 vs PC2
        plot(pca.EUR$PC1,
             pca.EUR$PC2,
             col = EUR.samples$Point_color,
             pch = EUR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC2', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
            # add legend
            legend('topright', legend = EUR.labels, col = EUR.colors, pch = EUR.pch, cex = 1, lty = 0, bty = 'n')
        
        # PC1 vs PC3
        plot(pca.EUR$PC1,
             pca.EUR$PC3,
             col = EUR.samples$Point_color,
             pch = EUR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC3', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC4
        plot(pca.EUR$PC1,
             pca.EUR$PC4,
             col = EUR.samples$Point_color,
             pch = EUR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC4', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC5
        plot(pca.EUR$PC1,
             pca.EUR$PC5,
             col = EUR.samples$Point_color,
             pch = EUR.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC5', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
    dev.off()
    
    # SAS
    pdf(paste(dir.output, "Plot_SNV_AND_SV_PCs_vs_Ancestry_SAS",".pdf", sep=""), width = 9, height = 9)
    par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses
    
        # PC1 vs PC2
        plot(pca.SAS$PC1,
             pca.SAS$PC2,
             col = SAS.samples$Point_color,
             pch = SAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC2', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
            # add legend
            legend('topright', legend = SAS.labels, col = SAS.colors, pch = SAS.pch, cex = 1, lty = 0, bty = 'n')
        
        # PC1 vs PC3
        plot(pca.SAS$PC1,
             pca.SAS$PC3,
             col = SAS.samples$Point_color,
             pch = SAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC3', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC4
        plot(pca.SAS$PC1,
             pca.SAS$PC4,
             col = SAS.samples$Point_color,
             pch = SAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC4', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
        # PC1 vs PC5
        plot(pca.SAS$PC1,
             pca.SAS$PC5,
             col = SAS.samples$Point_color,
             pch = SAS.samples$Point_pch,
             xlab = paste('PC1', sep=""),
             ylab = paste('PC5', sep=""),
             #xlim = c(-200, 200),
             #ylim = c(-90, 90),
             #main = "",
             cex.lab = pca.label.cex,
             cex.axis = pca.axis.cex
             )
        
    dev.off()
  
  

  

# Clean the environment
rm(list=ls())


