# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(BEDMatrix) # to load genotypes
library(ggplot2) # for violin plots
library (OmicCircos) # for circos plot

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/'




# SUMMARIZE L1 AND ALU SINGLETON COUNTS -------------------------------------------------------------------------------





# LOAD INPUT DATA

# Load sample meta-data
sample_info <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/integrated_call_samples_v3.20130502.ALL.panel.txt", header = T, row.names = 1, stringsAsFactors = F, sep = '\t')

    # Keep only the columns with data
    sample_info <- sample_info[, c('super_pop', 'pop', 'gender')] 
    
    # Update the column names
    colnames(sample_info) <- c('Super_Population', 'Population', 'Sex') 

    # Remove problematic samples from the meta-data
    samples_to_remove <- 'NA18498' # This sample has SV data, but no Phase3 hg38 SNV data
    sample_info <- sample_info[!(row.names(sample_info) %in% samples_to_remove), ]
    
# Specify the genotype BED files
plink_file_path.L1 <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Genotypes/Plink_L1_singletons/ALL_unique_LINE1.recode.plink.bed'
plink_file_path.Alu <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Genotypes/Plink_ALU_singletons/ALL_unique_ALU.recode.plink.bed'
plink_file_path.L1.TSD <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Genotypes/Plink_L1_singletons/With_TSD_ALL_unique_LINE1.recode.plink.bed'
plink_file_path.Alu.TSD <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Genotypes/Plink_ALU_singletons/With_TSD_ALL_unique_ALU.recode.plink.bed'

    # Make R object to stream snp genotypes into memory (rows are samples and columns are SVs). 
    binary_genotypes.L1 <- BEDMatrix(path = plink_file_path.L1, simple_names = T)
    binary_genotypes.Alu <- BEDMatrix(path = plink_file_path.Alu, simple_names = T)
    binary_genotypes.L1.TSD <- BEDMatrix(path = plink_file_path.L1.TSD, simple_names = T)
    binary_genotypes.Alu.TSD <- BEDMatrix(path = plink_file_path.Alu.TSD, simple_names = T)

    # Keep only the samples in the meta-data; also make sure sample order here matches the metadata table
    binary_genotypes.L1 <- binary_genotypes.L1[rownames(sample_info), ]  
    binary_genotypes.Alu <- binary_genotypes.Alu[rownames(sample_info), ]  
    binary_genotypes.L1.TSD <- binary_genotypes.L1.TSD[rownames(sample_info), ]  
    binary_genotypes.Alu.TSD <- binary_genotypes.Alu.TSD[rownames(sample_info), ]  
          
    # Check that only unique insertions are presents
    ncol(binary_genotypes.L1) == sum(colSums(binary_genotypes.L1))
    ncol(binary_genotypes.Alu) == sum(colSums(binary_genotypes.Alu))
    ncol(binary_genotypes.L1.TSD) == sum(colSums(binary_genotypes.L1.TSD))
    ncol(binary_genotypes.Alu.TSD) == sum(colSums(binary_genotypes.Alu.TSD))

    
    
    
    
# PREPARE SINGLETON SUMMARY DATA

# Sum the number of unique insertions in each sample (for L1, Alu, and L1+Alu)
sample_info$L1_unique_insertions <- rowSums(binary_genotypes.L1)
sample_info$Alu_unique_insertions <- rowSums(binary_genotypes.Alu)
sample_info$Combined_unique_insertions <- sample_info$L1_unique_insertions + sample_info$Alu_unique_insertions

sample_info$L1_TSD_unique_insertions <- rowSums(binary_genotypes.L1.TSD)
sample_info$Alu_TSD_unique_insertions <- rowSums(binary_genotypes.Alu.TSD)
sample_info$Combined_TSD_unique_insertions <- sample_info$L1_TSD_unique_insertions + sample_info$Alu_TSD_unique_insertions

# Add column to hold the unique insert info as a binary outcome in plink format (1 = no unique insertions and will be the control group in GWAS, 2 = at least 1 unique insertion and will be the case group in GWAS)
sample_info$L1_Plink_Binary <- 1
sample_info$Alu_Plink_Binary <- 1
sample_info$Combined_Plink_Binary <- 1

sample_info$L1_TSD_Plink_Binary <- 1
sample_info$Alu_TSD_Plink_Binary <- 1
sample_info$Combined_TSD_Plink_Binary <- 1

# Define the case samples (1 = no unique insertions and will be the control group in GWAS, 2 = at least 1 unique insertion and will be the case group in GWAS)
sample_info[sample_info$L1_unique_insertions > 0, 'L1_Plink_Binary'] <- 2
sample_info[sample_info$Alu_unique_insertions > 0, 'Alu_Plink_Binary'] <- 2
sample_info[sample_info$Combined_unique_insertions > 0, 'Combined_Plink_Binary'] <- 2

sample_info[sample_info$L1_TSD_unique_insertions > 0, 'L1_TSD_Plink_Binary'] <- 2
sample_info[sample_info$Alu_TSD_unique_insertions > 0, 'Alu_TSD_Plink_Binary'] <- 2
sample_info[sample_info$Combined_TSD_unique_insertions > 0, 'Combined_TSD_Plink_Binary'] <- 2

# Save the table
write.table(sample_info, file = paste(dir.output, "Singleton_Counts_Per_Sample_L1_and_Alu", '.txt', sep =""), sep = "\t" , row.names = T, col.names = NA, quote = F)







# GENERATE BARPLOTS FOR GWAS CASE/CONTROL FREQUENCIES

# To the sample info df, add a column with case/control status using the Combined L1+Alu counts
sample_info$GWAS_status <- 'controls'
sample_info[which(sample_info$Combined_unique_insertions > 0), 'GWAS_status'] <- 'cases'

# To the sample info df, add a column to use as a counter for cases and controls
sample_info$GWAS_counter <- 1

# Subset sample info by superpopulation
samples.AFR <- sample_info[which(sample_info$Super_Population == 'AFR'), ]
samples.AMR <- sample_info[which(sample_info$Super_Population == 'AMR'), ]
samples.EAS <- sample_info[which(sample_info$Super_Population == 'EAS'), ]
samples.EUR <- sample_info[which(sample_info$Super_Population == 'EUR'), ]
samples.SAS <- sample_info[which(sample_info$Super_Population == 'SAS'), ]

# Sum all the counter values by case/control status
GWAS_counts_AFR <- aggregate(samples.AFR$GWAS_counter, by = list(samples.AFR$GWAS_status), FUN = sum)
GWAS_counts_AMR <- aggregate(samples.AMR$GWAS_counter, by = list(samples.AMR$GWAS_status), FUN = sum)
GWAS_counts_EAS <- aggregate(samples.EAS$GWAS_counter, by = list(samples.EAS$GWAS_status), FUN = sum)
GWAS_counts_EUR <- aggregate(samples.EUR$GWAS_counter, by = list(samples.EUR$GWAS_status), FUN = sum)
GWAS_counts_SAS <- aggregate(samples.SAS$GWAS_counter, by = list(samples.SAS$GWAS_status), FUN = sum)
GWAS_counts_ALL <- aggregate(sample_info$GWAS_counter, by = list(sample_info$GWAS_status), FUN = sum)

# Generate bar plots
pdf(paste(dir.output, "Barplot_Case_Control_Counts_AFR", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(GWAS_counts_AFR$x ~ GWAS_counts_AFR$Group.1,
            ylab = "# of samples",
            ylim = c(0, 400),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # Add the exact counts
    text(0.70, 380, GWAS_counts_AFR[1, 2], cex = 1)
    text(1.90, 380, GWAS_counts_AFR[2, 2], cex = 1)
    
dev.off()


# Generate bar plots
pdf(paste(dir.output, "Barplot_Case_Control_Counts_AMR", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(GWAS_counts_AMR$x ~ GWAS_counts_AMR$Group.1,
            ylab = "# of samples",
            ylim = c(0, 400),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # Add the exact counts
    text(0.70, 380, GWAS_counts_AMR[1, 2], cex = 1)
    text(1.90, 380, GWAS_counts_AMR[2, 2], cex = 1)
    
dev.off()


# Generate bar plots
pdf(paste(dir.output, "Barplot_Case_Control_Counts_EAS", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(GWAS_counts_EAS$x ~ GWAS_counts_EAS$Group.1,
            ylab = "# of samples",
            ylim = c(0, 400),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # Add the exact counts
    text(0.70, 380, GWAS_counts_EAS[1, 2], cex = 1)
    text(1.90, 380, GWAS_counts_EAS[2, 2], cex = 1)
    
dev.off()


# Generate bar plots
pdf(paste(dir.output, "Barplot_Case_Control_Counts_EUR", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(GWAS_counts_EUR$x ~ GWAS_counts_EUR$Group.1,
            ylab = "# of samples",
            ylim = c(0, 400),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # Add the exact counts
    text(0.70, 380, GWAS_counts_EUR[1, 2], cex = 1)
    text(1.90, 380, GWAS_counts_EUR[2, 2], cex = 1)
    
dev.off()


# Generate bar plots
pdf(paste(dir.output, "Barplot_Case_Control_Counts_SAS", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(GWAS_counts_SAS$x ~ GWAS_counts_SAS$Group.1,
            ylab = "# of samples",
            ylim = c(0, 400),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # Add the exact counts
    text(0.70, 380, GWAS_counts_SAS[1, 2], cex = 1)
    text(1.90, 380, GWAS_counts_SAS[2, 2], cex = 1)
    
dev.off()

# Generate bar plots
pdf(paste(dir.output, "Barplot_Case_Control_Counts_ALL", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(GWAS_counts_ALL$x ~ GWAS_counts_ALL$Group.1,
            ylab = "# of samples",
            ylim = c(0, 1600),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # Add the exact counts
    text(0.70, 1500, GWAS_counts_ALL[1, 2], cex = 1)
    text(1.90, 1500, GWAS_counts_ALL[2, 2], cex = 1)
    
dev.off()







# GENERATE VIOLIN PLOTS FOR SINGLETON DISTRIBUTION AMONG CASES

# Subset the cases in each superpopulation
cases.ALL <- sample_info[which(sample_info$GWAS_status == 'cases'), ]

# Order factors levels
cases.ALL$Super_Population <- factor(cases.ALL$Super_Population, levels = c("AFR", "EAS", "EUR", "SAS", "AMR"))

# Generate violin plots
pdf(paste(dir.output, "Violinplot_Singleton_Distribution_in_cases", ".pdf", sep=""), width = 6, height = 6)

ggplot(cases.ALL, aes(x = Super_Population, y = Combined_unique_insertions)) + 
        geom_violin(aes(fill = Super_Population), trim = FALSE) +
        geom_boxplot(width = 0.1) +
        scale_y_continuous(breaks = seq(0, 20, by = 2)) +
        scale_fill_manual(values=c(AFR = "#4F73FF", EAS = "#CA2A64", EUR ="#FC4C1A", SAS = "#FF941A", AMR = "#5F4BE7")) +
        theme_classic() +
        theme(legend.position = "none")

dev.off()





# VISUALIZE SINGLETON FREQUENCIES AND LOCATIONS -------------------------------------------------------------------------------





# PREPARE SINGLETON LOCATIONS FOR CIRCOS

# Load singleton data from plink BIM files
L1_singleton_info <- read.csv('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Genotypes/Plink_L1_singletons/ALL_unique_LINE1.recode.plink.bim', header = F, stringsAsFactors = F, sep = '\t')
Alu_singleton_info <- read.csv('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Genotypes/Plink_ALU_singletons/ALL_unique_ALU.recode.plink.bim', header = F, stringsAsFactors = F, sep = '\t')

# Rename columns.
colnames(L1_singleton_info) <- c('chr', 'SV_ID', 'Morgan_Pos', 'pos', 'ALT', 'REF')
colnames(Alu_singleton_info) <- c('chr', 'SV_ID', 'Morgan_Pos', 'pos', 'ALT', 'REF')

# Make function to add each sample and superpopulation corresponding to each unique insertion
Map_singleton_to_sample <- function(Singleton_info, Singleton_genotypes, sample_metadata){
  
  # Add a column to hold the sample ID
  Singleton_info$Sample <- c('')
  
  # Run for loop over each SV
  for (SV_index in 1:nrow(Singleton_info)) {
    
    # Extract the name of the ith SV
    ith_SV_name <- Singleton_info[SV_index, 'SV_ID']

    # Find the sample that has that SV, and add the sample name to the SV info table
    Singleton_info[SV_index, 'Sample'] <- names(which(Singleton_genotypes[, ith_SV_name] == 1))
    
  } # END for loop
  
  # Map each sample to a superpopulation using a metadata table, and add that info to the SV info table
  Singleton_info$Super_Population <- sample_metadata[Singleton_info$Sample, 'Super_Population']
  
  # Output the updated SV info table
  return(Singleton_info)
     
}

# Run function to add each sample and superpopulation corresponding to each unique insertion
L1_singleton_info <- Map_singleton_to_sample(Singleton_info = L1_singleton_info, Singleton_genotypes = binary_genotypes.L1, sample_metadata = sample_info)
Alu_singleton_info <- Map_singleton_to_sample(Singleton_info = Alu_singleton_info, Singleton_genotypes = binary_genotypes.Alu, sample_metadata = sample_info)

# Save the tables before modifying further
write.table(L1_singleton_info, file = paste(dir.output, "Singleton_Details_L1_only", '.txt', sep =""), sep = "\t" , row.names = F, col.names = TRUE, quote = F)
write.table(Alu_singleton_info, file = paste(dir.output, "Singleton_Details_Alu_only", '.txt', sep =""), sep = "\t" , row.names = F, col.names = TRUE, quote = F)
write.table(rbind(L1_singleton_info, Alu_singleton_info), file = paste(dir.output, "Singleton_Details_L1_and_Alu_Combined", '.txt', sep =""), sep = "\t" , row.names = F, col.names = TRUE, quote = F)

# Append the 'chr' label to chromosomes
L1_singleton_info$chr <- paste('chr', L1_singleton_info$chr, sep = '')
Alu_singleton_info$chr <- paste('chr', Alu_singleton_info$chr, sep = '')
  
# Add a column with a numeric value that decides the length of each peak in the circos plot (this will be arbitrarily set to 1)
L1_singleton_info$Intensity <- 1
Alu_singleton_info$Intensity <- 1

# Re-arrange columns to match circos input
circos.L1 <- L1_singleton_info[, c('chr', 'pos', 'Super_Population', 'Intensity')]
circos.Alu <- Alu_singleton_info[, c('chr', 'pos', 'Super_Population', 'Intensity')]

# Update the colnames
colnames(circos.L1) <- c('chromosome_name', 'start_position', 'Super_Population', 'Intensity')
colnames(circos.Alu) <- c('chromosome_name', 'start_position', 'Super_Population', 'Intensity')

# Generate a table with L1 and Alu singletons combined
circos.combined <- rbind(circos.L1, circos.Alu)
  
    
            




# PREPARE CIRCOS BACKBONE USING CYTOGENETIC DATA
    
# Load cytogenetic data
hg38_cytobands <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/cytoBand_hg38.txt", header = FALSE, stringsAsFactors = F, sep = '\t', row.names = NULL)

# Assign column names
colnames(hg38_cytobands) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain')

# Change gieStain to colors
hg38_cytobands[hg38_cytobands$gieStain == c('gneg'), 'gieStain'] <- c('gray88')
hg38_cytobands[hg38_cytobands$gieStain == c('gpos25'), 'gieStain'] <- c('gray68')
hg38_cytobands[hg38_cytobands$gieStain == c('gpos50'), 'gieStain'] <- c('gray58')
hg38_cytobands[hg38_cytobands$gieStain == c('gpos75'), 'gieStain'] <- c('gray48')
hg38_cytobands[hg38_cytobands$gieStain == c('gpos100'), 'gieStain'] <- c('gray38')
hg38_cytobands[hg38_cytobands$gieStain == c('acen'), 'gieStain'] <- c('red') # "acen" is centromeric
hg38_cytobands[hg38_cytobands$gieStain == c('gvar'), 'gieStain'] <- c('gray88') # "gvar" bands tend to be heterochomatin, either pericentric or telomeric
hg38_cytobands[hg38_cytobands$gieStain == c('stalk'), 'gieStain'] <- c('gray88') # "stalk" refers to the short arm of acrocentric chromosomes chr13,14,15,21,22

# Make a vector of the chromosomes to keep
chr_to_keep <- paste(c('chr'), 1:22, sep = '')

# Define indices of the cytogenetic input table where the chromosomes labels match the chromosomes to keep
chr_to_keep_indices <- hg38_cytobands$chrom %in% chr_to_keep

# Extract the desired info from the cytogenetic band table
hg38_cytobands <- hg38_cytobands[chr_to_keep_indices, ]

# Convert backbone linear coordinates to radial coordinates
chr_radial_coord <- segAnglePo(seg.dat = hg38_cytobands, seg = chr_to_keep, angle.start = 0, angle.end = 350)
            
            
            
            


            
# GENERATE THE CIRCOS PLOTS

# L1 singletons circos plot
pdf(paste(dir.output, "Circos_Plot_Singletons_L1_only", ".pdf", sep=""))
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab = "")

    # Plot the backbone (no color, but necessary to keep the chr names)
    circos(cir = chr_radial_coord, type = 'chr', col = FALSE, print.chr.lab = TRUE, R = 350, W = 100, scale = FALSE)
    
    # Overlay the cytobands over the backbone
    circos(mapping = hg38_cytobands, cir = chr_radial_coord, type = 'arc2', col.v = 2, col = hg38_cytobands$gieStain, print.chr.lab = TRUE, R = 350, W = 0, lwd = 10, cex = 1, scale = FALSE)
    
    # Overlay the EAS singletons
    circos(mapping = circos.L1[which(circos.L1$Super_Population == 'EAS'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#DC267F"), 
           R = 300, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the AFR singletons
    circos(mapping = circos.L1[which(circos.L1$Super_Population == 'AFR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#648FFF"), 
           R = 250, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the SAS singletons
    circos(mapping = circos.L1[which(circos.L1$Super_Population == 'SAS'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#FFB000"), 
           R = 200, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the EUR singletons
    circos(mapping = circos.L1[which(circos.L1$Super_Population == 'EUR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#FE6100"), 
           R = 150, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the AMR singletons
    circos(mapping = circos.L1[which(circos.L1$Super_Population == 'AMR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#785EF0"), 
           R = 100, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
# End plot
dev.off()   
            
            

    

# Alu singletons circos plot
pdf(paste(dir.output, "Circos_Plot_Singletons_Alu_only", ".pdf", sep=""))
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab = "")

    # Plot the backbone (no color, but necessary to keep the chr names)
    circos(cir = chr_radial_coord, type = 'chr', col = FALSE, print.chr.lab = TRUE, R = 350, W = 100, scale = FALSE)
    
    # Overlay the cytobands over the backbone
    circos(mapping = hg38_cytobands, cir = chr_radial_coord, type = 'arc2', col.v = 2, col = hg38_cytobands$gieStain, print.chr.lab = TRUE, R = 350, W = 0, lwd = 10, cex = 1, scale = FALSE)
    
    # Overlay the AFR singletons
    circos(mapping = circos.Alu[which(circos.Alu$Super_Population == 'AFR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#648FFF"), 
           R = 300, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the EAS singletons
    circos(mapping = circos.Alu[which(circos.Alu$Super_Population == 'EAS'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#DC267F"), 
           R = 250, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the EUR singletons
    circos(mapping = circos.Alu[which(circos.Alu$Super_Population == 'EUR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#FE6100"), 
           R = 200, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the SAS singletons
    circos(mapping = circos.Alu[which(circos.Alu$Super_Population == 'SAS'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#FFB000"), 
           R = 150, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the AMR singletons
    circos(mapping = circos.Alu[which(circos.Alu$Super_Population == 'AMR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#785EF0"), 
           R = 100, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
# End plot
dev.off()       
            
            

    

# Combined singletons circos plot
pdf(paste(dir.output, "Circos_Plot_Singletons_L1_and_Alu_Combined", ".pdf", sep=""))
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type = "n", axes = FALSE, xlab = "", ylab = "")

    # Plot the backbone (no color, but necessary to keep the chr names)
    circos(cir = chr_radial_coord, type = 'chr', col = FALSE, print.chr.lab = TRUE, R = 350, W = 100, scale = FALSE)
    
    # Overlay the cytobands over the backbone
    circos(mapping = hg38_cytobands, cir = chr_radial_coord, type = 'arc2', col.v = 2, col = hg38_cytobands$gieStain, print.chr.lab = TRUE, R = 350, W = 0, lwd = 10, cex = 1, scale = FALSE)
    
    # Overlay the AFR singletons
    circos(mapping = circos.combined[which(circos.combined$Super_Population == 'AFR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#648FFF"), 
           R = 300, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the EAS singletons
    circos(mapping = circos.combined[which(circos.combined$Super_Population == 'EAS'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#DC267F"), 
           R = 250, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the EUR singletons
    circos(mapping = circos.combined[which(circos.combined$Super_Population == 'EUR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#FE6100"), 
           R = 200, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the SAS singletons
    circos(mapping = circos.combined[which(circos.combined$Super_Population == 'SAS'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#FFB000"), 
           R = 150, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
    # Overlay the AMR singletons
    circos(mapping = circos.combined[which(circos.combined$Super_Population == 'AMR'), ], 
           cir = chr_radial_coord, 
           type = 'b3', 
           col = c("#785EF0"), 
           R = 100, 
           W = 50, 
           lwd = 1, 
           cex = 1, 
           cutoff = 0, 
           B = TRUE, 
           print.chr.lab = FALSE, 
           scale = FALSE)
    
# End plot
dev.off()       







# GENERATE BARPLOTS FOR SINGLETON FREQUENCY IN EACH SUPERPOPULATION

# Aggregate singleton counts by superpopulation
L1_aggregate_counts <- aggregate(circos.L1$Intensity, by = list(circos.L1$Super_Population), FUN = sum)
Alu_aggregate_counts <- aggregate(circos.Alu$Intensity, by = list(circos.Alu$Super_Population), FUN = sum)
Combined_aggregate_counts <- aggregate(circos.combined$Intensity, by = list(circos.combined$Super_Population), FUN = sum)

# Update column names
colnames(L1_aggregate_counts) <- c('Super_Population', 'Counts')
colnames(Alu_aggregate_counts) <- c('Super_Population', 'Counts')
colnames(Combined_aggregate_counts) <- c('Super_Population', 'Counts')

# Order results from largest to smallest counts
L1_aggregate_counts <- L1_aggregate_counts[order(L1_aggregate_counts$Counts, decreasing = TRUE), ]
Alu_aggregate_counts <- Alu_aggregate_counts[order(Alu_aggregate_counts$Counts, decreasing = TRUE), ]
Combined_aggregate_counts <- Combined_aggregate_counts[order(Combined_aggregate_counts$Counts, decreasing = TRUE), ]

# Define group levels (needed to specify order in the plot)
L1_aggregate_counts$Super_Population <- factor(L1_aggregate_counts$Super_Population, levels = L1_aggregate_counts$Super_Population)
Alu_aggregate_counts$Super_Population <- factor(Alu_aggregate_counts$Super_Population, levels = Alu_aggregate_counts$Super_Population)
Combined_aggregate_counts$Super_Population <- factor(Combined_aggregate_counts$Super_Population, levels = Combined_aggregate_counts$Super_Population)

# Generate L1 singleton bar plots
pdf(paste(dir.output, "Barplot_Singleton_Counts_Per_Superpopulation_L1_only", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(L1_aggregate_counts$Counts ~ L1_aggregate_counts$Super_Population,
            ylab = "# of L1 Singletons",
            ylim = c(0, 400),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
dev.off()
                       
# Generate Alu singleton bar plots
pdf(paste(dir.output, "Barplot_Singleton_Counts_Per_Superpopulation_Alu_only", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(Alu_aggregate_counts$Counts ~ Alu_aggregate_counts$Super_Population,
            ylab = "# of Alu Singletons",
            ylim = c(0, 700),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
dev.off()      
              
# Generate Combined singleton bar plots
pdf(paste(dir.output, "Barplot_Singleton_Counts_Per_Superpopulation_L1_and_Alu_Combined", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(Combined_aggregate_counts$Counts ~ Combined_aggregate_counts$Super_Population,
            ylab = "# of Alu and L1 Singletons",
            ylim = c(0, 1000),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
dev.off()               
    

   



       
# SESSION INFO -------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_Singleton_Frequencies_and_Distributions.txt", sep =""))
sessionInfo()
sink()      



# Clean the environment
rm(list=ls())



