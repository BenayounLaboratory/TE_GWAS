# Set strings as factors
options(stringsAsFactors = F)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/1_Population_Structure_Analysis_SNV_and_SV/Population_Structure_Analysis/'



# Load GEUVADIS metadata
GEUV_sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/E-GEUV-1-unique_and_in_Phase3_snps_TE_SVs.txt", header = T, stringsAsFactors = F, sep = '\t')
rownames(GEUV_sample_info) <- GEUV_sample_info$Source.Name

# Extract EUR and YRI metadata
EUR.samples <- GEUV_sample_info[which(GEUV_sample_info$Characteristics.ancestry.category. != 'Yoruba'), ]
YRI.samples <- GEUV_sample_info[which(GEUV_sample_info$Characteristics.ancestry.category. == 'Yoruba'), ]

# Load PCA eigenvectors.
pca.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/1_Population_Structure_Analysis_SNV_and_SV/Population_Structure_Genotypes/Plink_EUR_PCA/Final_SNV_SV_EUR_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")
pca.YRI <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/1_Population_Structure_Analysis_SNV_and_SV/Population_Structure_Genotypes/Plink_YRI_PCA/Final_SNV_SV_YRI_PCA.eigenvec", header = FALSE, row.names = 1, sep = "")

    # Match the order of the sample names
    pca.EUR <- pca.EUR[rownames(EUR.samples), ]
    pca.YRI <- pca.YRI[rownames(YRI.samples), ]
    
    # Clean PCA df: remove sample names column and assign PC labels
    pca.EUR <- pca.EUR[, -c(1)]
    pca.YRI <- pca.YRI[, -c(1)]
    
    colnames(pca.EUR) <- paste('PC', 1:ncol(pca.EUR), sep = '')
    colnames(pca.YRI) <- paste('PC', 1:ncol(pca.YRI), sep = '')
    
# Load PCA eigenvalues.
eigenval.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/1_Population_Structure_Analysis_SNV_and_SV/Population_Structure_Genotypes/Plink_EUR_PCA/Final_SNV_SV_EUR_PCA.eigenval", header = FALSE)
eigenval.YRI <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/1_Population_Structure_Analysis_SNV_and_SV/Population_Structure_Genotypes/Plink_YRI_PCA/Final_SNV_SV_YRI_PCA.eigenval", header = FALSE)

    # Calculate percent variance explained (pve)
    pve.EUR <- 100 * (eigenval.EUR/sum(eigenval.EUR)) 
    pve.YRI <- 100 * (eigenval.YRI/sum(eigenval.YRI)) 

# Add a 'color' column to the sample tables
EUR.samples$Point_color <- '#FE6100'
YRI.samples$Point_color <- '#5A81E6'

    # Update the point colors
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'British'), 'Point_color'] <- '#190A00'
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Utah'), 'Point_color'] <- '#FE9859'
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Finnish'), 'Point_color'] <- '#7F3100'
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Tuscan'), 'Point_color'] <- '#FFD8BF'

# Add a 'shape' column to the sample table
EUR.samples$Point_pch <- as.numeric('16')
YRI.samples$Point_pch <- as.numeric('3')

    # Update the point shapes
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'British'), 'Point_pch'] <- 15
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Utah'), 'Point_pch'] <- 18
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Finnish'), 'Point_pch'] <- 16
    EUR.samples[which(EUR.samples$Characteristics.ancestry.category. == 'Tuscan'), 'Point_pch'] <- 3

# Define plot and legend parameters
EUR.labels <- c('GBR', 'CEU', 'FIN', 'TSI')
EUR.colors <- c('#190A00', '#FE9859', '#7F3100',  '#FFD8BF')
EUR.pch <- c(15, 18, 16, 3)

YRI.labels <- c('YRI')
YRI.colors <- c('#5A81E6')
YRI.pch <- c(3)

pca.label.cex <- 1
pca.axis.cex <- 1





# Plot PVE
pdf(paste(dir.output, "Plot_Percent_Variance_Explained_SNVs_and_SVs", ".pdf", sep=""))
par(mfrow=c(2,1))
    barplot(t(pve.EUR), names.arg = rownames(pve.EUR), xlab = "Principal component", ylab = 'Percent variance explained', main = 'EUR Genotype PCA', las = 2)
    barplot(t(pve.YRI), names.arg = rownames(pve.YRI), xlab = "Principal component", ylab = 'Percent variance explained', main = 'YRI Genotype PCA', las = 2)
dev.off()

# Plot PCs vs Ancestry
pdf(paste(dir.output, "Plot_Principal_Components_vs_Ancestry_SNVs_and_SVs",".pdf", sep=""), width = 7, height = 7)
par(mfrow=c(2,2)) # Added to produce images of similar shape across analyses

    # EUR: PC1 vs PC2
    plot(pca.EUR$PC1,
         pca.EUR$PC2,
         col = EUR.samples$Point_color,
         pch = EUR.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC2', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "EUR Population"
         )
    
        # add legend
        legend('topright', legend = EUR.labels, col = EUR.colors, pch = EUR.pch, cex = 1, lty = 0, bty = 'n')
    
    # EUR: PC1 vs PC3
    plot(pca.EUR$PC1,
         pca.EUR$PC3,
         col = EUR.samples$Point_color,
         pch = EUR.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC3', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "EUR Population"
         )
    
    # YRI: PC1 vs PC2
    plot(pca.YRI$PC1,
         pca.YRI$PC2,
         col = YRI.samples$Point_color,
         pch = YRI.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC2', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "YRI Population"
         )
        
        # add a legend
        legend('topright', legend = YRI.labels, col = YRI.colors, pch = YRI.pch, cex = 1, lty = 0, bty = 'n')
    
    # YRI: PC1 vs PC3
    plot(pca.YRI$PC1,
         pca.YRI$PC3,
         col = YRI.samples$Point_color,
         pch = YRI.samples$Point_pch,
         xlab = paste('PC1', sep=""),
         ylab = paste('PC3', sep=""),
         #xlim = c(-200, 200),
         #ylim = c(-90, 90),
         cex.lab = pca.label.cex,
         cex.axis = pca.axis.cex,
         main = "YRI Population"
         )
    
    
dev.off()

    
  

  
# Combine EUR/YRI PCs and save
pca.EUR.YRI <- rbind(pca.EUR, pca.YRI)
write.table(pca.EUR.YRI, file = paste(dir.output, "COVARIATE_Population_Structure_PCA_SNVs_and_SVs", ".txt", sep =""), sep = "\t" , row.names = T, col.names = NA, quote=F)
  
  
  


# Clean the environment
rm(list=ls())

