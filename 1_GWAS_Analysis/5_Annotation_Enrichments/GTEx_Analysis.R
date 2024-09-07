# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(pheatmap) # for fread and fsave

# Specify directory for annotated SNVs
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/GTEx_Expression_Analysis/'





# LOAD AND PROCESS DATA

# Load genelist for SNVs
genes.greenlist.SNVs <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_Significant_Variants_Greenlist_SNVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)

# Load GTEx median TPM per tissue expression file
median.GTEx <- read.delim(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/GTEx_Analysis_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', header = TRUE, sep = "\t", skip = 2)

    # Only keep GTEx expression values for genes associated to GWAS greenlist SNVs
    median.GTEx <- median.GTEx[which(median.GTEx$Description %in% genes.greenlist.SNVs$V1), ]
    
    # Check for gene symbol duplicates
    sum(duplicated(median.GTEx$Description)) # 1 gene symbol duplicated.
    
    # Remove gene symbol duplicates from GTEx expression
    median.GTEx <- median.GTEx[!duplicated(median.GTEx$Description), ]
    
    # Assign gene symbols to rows
    rownames(median.GTEx) <- median.GTEx$Description
    
    # Remove gene name columns
    median.GTEx <- median.GTEx[, -c(1,2)]
    
    # Remove rows with all zeros
    median.GTEx <- median.GTEx[rowSums(median.GTEx) > 0,] # 8 genes removed





# GENERATE HEATMAPS
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)
    
# Generate preliminary heatmap
p1 <- pheatmap(median.GTEx, 
               cluster_rows = TRUE, 
               show_rownames = FALSE,
               fontsize_row = 1,
               cluster_cols = TRUE, 
               show_colnames = FALSE,
               scale = 'row',
               border_color = NA,
               #annotation_col = annot_column,
               #annotation_colors = annot_color,
               color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
               breaks = breaksList,
               main = paste('', sep = ''))

dev.off()

# Extract pheatmap column order
my.col.order <- p1$tree_col$order
my.col.order <- colnames(median.GTEx[, my.col.order])
  
# Extract pheatmap row order
my.row.order <- p1$tree_row$order
my.row.order <- rownames(median.GTEx[my.row.order, ])  
    
# Update the gene expression matrix order
median.GTEx <- median.GTEx[my.row.order, my.col.order]

    # Find indices for gonad and brain tissues
    indices.gonad <- which(colnames(median.GTEx) %in% c('Testis', 'Ovary'))
    indices.brain <- which(colnames(median.GTEx) %in% c('Brain...Cerebellar.Hemisphere', 'Brain...Cerebellum', 'Brain...Frontal.Cortex..BA9.', 'Brain...Anterior.cingulate.cortex..BA24.', 'Brain...Cortex', 'Brain...Spinal.cord..cervical.c.1.', 'Brain...Nucleus.accumbens..basal.ganglia.', 'Brain...Caudate..basal.ganglia.', 'Brain...Putamen..basal.ganglia.', 'Brain...Hypothalamus', 'Brain...Substantia.nigra', 'Brain...Amygdala', 'Brain...Hippocampus'))
    indices.combined <- c(indices.gonad, indices.brain)
    
    # Reorder expression matrix so gonads and brain are first (this will simplify separating them during figure creation)
    median.GTEx <- cbind(median.GTEx[, indices.combined], median.GTEx[, -c(indices.combined)])
  

# Save updated table
write.table(median.GTEx, file = paste(dir.output, "GTEx_Greenlist_SNV_Associated_Gene_Expression", '.txt', sep=""), sep = "\t" , row.names = T, col.names = NA, quote = F)         

# Save heatmap
pdf(paste(dir.output, "Plot_Heatmap_Greenlist_SNV_Associated_Gene_Expression", ".pdf", sep=""), height = 8, width = 10)

  pheatmap(median.GTEx, 
               cluster_rows = FALSE, 
               show_rownames = FALSE,
               fontsize_row = 1,
               cluster_cols = FALSE, 
               show_colnames = FALSE,
               scale = 'row',
               border_color = NA,
               #annotation_col = annot_column,
               #annotation_colors = annot_color,
               color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
               breaks = breaksList,
               main = paste('', sep = ''))
dev.off()
dev.off()



# Clean the environment
rm(list=ls())       
        

