# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(data.table) # for fread and fsave
library(splitstackshape) # to split concatenated genes/TEs

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Calculate_Enrichments_SNVs_functions.R") 

# Specify directory for annotated SNVs
dir.output.green.SNV <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Enrichment_Results_SNVs/Greenlist_SNVs/'





# LOAD SNV DATA

# Load all/background annotated GWAS SNP lists
all.snvs <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/All_Annotated_GWAS_Background_SNVs_only.csv', header = TRUE, sep = ",", data.table = FALSE)

    # Check whether there are duplicate RsIDs
    sum(duplicated(all.snvs$RsID)) # 0
    
# Load significant annotated GWAS SNV/SV list
sig.variants.all <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/All_Annotated_GWAS_Significant_SNVs_AND_SVs.csv', header = TRUE, sep = ",")

    # ONLY keep the SNVs
    sig.snvs.all <- sig.variants.all[which(sig.variants.all$SNP %in% all.snvs$RsID), ]
    
    # Check whether there are duplicate RsIDs
    sum(duplicated(sig.snvs.all$RsID)) # 0
    
# Make a second sig snp list with greenlist SNVs only
    
    # make a second df to hold non-blacklisted SNVs
    sig.snvs.greenlist <- sig.snvs.all
    
    # Remove blacklisted SNPs
    sig.snvs.greenlist <- sig.snvs.greenlist[which(is.na(sig.snvs.greenlist$Blacklist_Label)), ]
    
# Make a df to hold enrichment pvalues and FDR
all.enrichments <- data.frame(p = rep(NA, 18), FDR = rep(NA, 18))
    
    # Assign rownames
    rownames(all.enrichments) <- c('Blacklist', 'SV_Hotspot', 'Segmental_Duplication', 'Transposition_Regulators', 'Expression_Regulators', 'Methyltransferase', 'RNA_Modification', 'Bravo_Regulators', 'AluJ', 'AluS', 'AluY', 'L1M', 'L1P', 'L1PA', 'Intact-Full_L1', 'Nonintact-Full_L1', 'Intact_ORF2_L1', 'ENCODE_cCREs')
    
    
    
    
    
# CHECK ENRICHMENT OF ENCODE BLACKLIST ANNOTATIONS
    
# Calculate proportions in the significant SNP list
    
    # number of snps with 'high signal region' or 'low mapability' (i.e. do not have the NA label)
    sig.in.blacklist <- length(which(!is.na(sig.snvs.all$Blacklist_Label)))
    
    # number of snps without 'high signal region' or 'low mapability'
    sig.notin.blacklist <- nrow(sig.snvs.all) - sig.in.blacklist
    
    # fraction of snps with 'high signal region
    sig.fraction.blacklist <- sig.in.blacklist / (sig.in.blacklist + sig.notin.blacklist)
    
# Calculate proportions in the background SNP list
    
    # number of snps with 'high signal region' or 'low mapability' (i.e. do not have the NA label)
    background.in.blacklist <- length(which(!is.na(all.snvs$Blacklist_Label)))
    
    # number of snps without 'high signal region' or 'low mapability'
    background.notin.blacklist <- nrow(all.snvs) - background.in.blacklist
    
    # fraction of snps with 'high signal region
    background.fraction.blacklist <- background.in.blacklist / (background.in.blacklist + background.notin.blacklist)
    
# Collect fractions in a df to use for plot
fraction.blacklisted <- data.frame(Fraction = c(sig.fraction.blacklist, background.fraction.blacklist), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.blacklist <- as.matrix(data.frame(Not_Blacklisted = c(sig.notin.blacklist, background.notin.blacklist), Blacklisted = c(sig.in.blacklist, background.in.blacklist)))

# Run Fischer's Exact Test
test.blacklist <- fisher.test(my.contingency.blacklist)
  
# Add pvalue to the table
all.enrichments['Blacklist', 'p'] <- test.blacklist$p.value  

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_Blacklisted_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.blacklisted$Fraction ~ fraction.blacklisted$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs in ENCODE Blacklist v2',
            ylim = c(0, 1),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, .88, paste('p = ', formatC(test.blacklist$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.85,
             x1 = 1.88, y1 = 0.85)
    
dev.off()

    



# CHECK ENRICHMENT OF SV HOTSPOT (Lin 2019) ANNOTATIONS
  
# Calculate fractions with annotation and Fischer test results
analysis.SV_Hotspot <- Annotation_enrichment_sig_vs_background(annotated_background_variants = all.snvs, 
                                                                   annotated_significant_variants = sig.snvs.greenlist, 
                                                                   annotation_column_name = 'SV_Hotspot_Label', 
                                                                   annotation_label = 'Lin_2019')
    
# Extract fractions and Fischer test results
fraction.SV_hotspot <- analysis.SV_Hotspot[[1]]
test.SV_hotspot <- analysis.SV_Hotspot[[2]]
  
# Add pvalue to the table
all.enrichments['SV_Hotspot', 'p'] <- test.SV_hotspot$p.value  

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_SV-Hotspot_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.SV_hotspot$Fraction ~ fraction.SV_hotspot$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs in SV Hotspots (Lin 2019)',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.SV_hotspot$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()

    



# CHECK ENRICHMENT OF SEGMENTAL DUPLICATIONS (BAILEY 2001) ANNOTATIONS
    
# Calculate fractions with annotation and Fischer test results
analysis.SegmentalDup <- Annotation_enrichment_sig_vs_background(annotated_background_variants = all.snvs, 
                                                                   annotated_significant_variants = sig.snvs.greenlist, 
                                                                   annotation_column_name = 'Segmental_Duplication', 
                                                                   annotation_label = 'Bailey_2001')
    
# Extract fractions and Fischer test results
fractions.SegmentalDup <- analysis.SegmentalDup[[1]]
test.SegmentalDup <- analysis.SegmentalDup[[2]]
  
# Add pvalue to the table
all.enrichments['Segmental_Duplication', 'p'] <- test.SegmentalDup$p.value  

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_Segmental_Duplication_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fractions.SegmentalDup$Fraction ~ fractions.SegmentalDup$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs in Segmental Duplication (Bailey 2001)',
            ylim = c(0, 0.40),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.39, paste('p = ', formatC(test.SegmentalDup$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.38,
             x1 = 1.88, y1 = 0.38)
    
dev.off()





# CHECK ENRICHMENT OF KNOWN L1 TRANSPOSITION REGULATORS

# Load file with known L1 transposition regulators. Use the gene symbols as the reference list.
known.transposition <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Liu_2018_Transposition_Regulators/L1_regulators_Primary_Screen.txt', header = TRUE, sep = "\t")
known.transposition <- known.transposition$Symbol

# Calculate fractions with annotation and Fischer test results
analysis.transposition <- Gene_overlap_enrichment_sig_vs_background(annotated_background_variants = all.snvs[, c('RsID', 'Nearby_Genes')], 
                                                                    annotated_significant_variants = sig.snvs.greenlist[, c('SNP', 'Nearby_Genes')], 
                                                                    genes_of_interest = known.transposition)

# Extract fractions, Fischer test results, and variants overlapping genes of interest
fraction.transposition <- analysis.transposition[[1]]
test.transposition <- analysis.transposition[[2]]
sig.variants.transposition <- analysis.transposition[[3]]
    
# Add pvalue to the table
all.enrichments['Transposition_Regulators', 'p'] <- test.transposition$p.value

# Save table of variants overlapping genes of interest
write.table(sig.variants.transposition, file = paste(dir.output.green.SNV, "Variants_overlapping_L1_Transposition_Regulators", ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA", quote = F)

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_Transposition_Regulators_Near_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.transposition$Fraction ~ fraction.transposition$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs near L1 transposition regulators (Liu et al. 2018)',
            ylim = c(0, 0.2),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, .19, paste('p = ', formatC(test.transposition$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF KNOWN L1 EXPRESSION REGULATORS

# Load file with known L1 transposition regulators. Use the gene symbols as the reference list.
known.expression <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Li_2024_L1_Expression_Regulators/Li_2024_L1_Expression_Regulators_Both_Activators_Suppressors.txt', header = TRUE, sep = "\t")
known.expression <- known.expression$id

# Calculate fractions with annotation and Fischer test results
analysis.expression <- Gene_overlap_enrichment_sig_vs_background(annotated_background_variants = all.snvs[, c('RsID', 'Nearby_Genes')], 
                                                                    annotated_significant_variants = sig.snvs.greenlist[, c('SNP', 'Nearby_Genes')], 
                                                                    genes_of_interest = known.expression)

# Extract fractions, Fischer test results, and variants overlapping genes of interest
fraction.expression <- analysis.expression[[1]]
test.expression <- analysis.expression[[2]]
sig.variants.expression <- analysis.expression[[3]]
    
# Add pvalue to the table
all.enrichments['Expression_Regulators', 'p'] <- test.expression$p.value

# Save table of variants overlapping genes of interest
write.table(sig.variants.expression, file = paste(dir.output.green.SNV, "Variants_overlapping_L1_Expression_Regulators", ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA", quote = F)

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_Expression_Regulators_Near_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.expression$Fraction ~ fraction.expression$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs near L1 expression regulators (Li et al. 2024)',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, .19, paste('p = ', formatC(test.expression$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF HISTONE METHYLTRANSFERASE GENES

# Load file with known genes. 
known.methyltransferase <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Other_Potential_TE_Regulators/MSigDB_GOMF_HISTONE_METHYLTRANSFERASE_ACTIVITY.v2023.2.Hs.txt', header = FALSE, sep = "\t")
known.methyltransferase <- known.methyltransferase$V1

# Calculate fractions with annotation and Fischer test results
analysis.methyl <- Gene_overlap_enrichment_sig_vs_background(annotated_background_variants = all.snvs[, c('RsID', 'Nearby_Genes')], 
                                                                    annotated_significant_variants = sig.snvs.greenlist[, c('SNP', 'Nearby_Genes')], 
                                                                    genes_of_interest = known.methyltransferase)

# Extract fractions, Fischer test results, and variants overlapping genes of interest
fraction.methyl <- analysis.methyl[[1]]
test.methyl <- analysis.methyl[[2]]
sig.variants.methyl <- analysis.methyl[[3]]
    
# Add pvalue to the table
all.enrichments['Methyltransferase', 'p'] <- test.methyl$p.value

# Save table of variants overlapping genes of interest
write.table(sig.variants.methyl, file = paste(dir.output.green.SNV, "Variants_overlapping_Histone_Methyltransferases", ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA", quote = F)

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_Histone_Methyltransferases_Near_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.methyl$Fraction ~ fraction.methyl$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs near histone methyltransferases (GO:0042054)',
            ylim = c(0, 0.10),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, .09, paste('p = ', formatC(test.methyl$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.08,
             x1 = 1.88, y1 = 0.08)
    
dev.off()





# CHECK ENRICHMENT OF RNA MODIFICATION GENES

# Load file with known genes. 
known.RNA_modification <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Other_Potential_TE_Regulators/MSigDB_GOBP_RNA_MODIFICATION.v2023.2.Hs.txt', header = FALSE, sep = "\t")
known.RNA_modification <- known.RNA_modification$V1

# Calculate fractions with annotation and Fischer test results
analysis.RNA_modification <- Gene_overlap_enrichment_sig_vs_background(annotated_background_variants = all.snvs[, c('RsID', 'Nearby_Genes')], 
                                                                       annotated_significant_variants = sig.snvs.greenlist[, c('SNP', 'Nearby_Genes')], 
                                                                       genes_of_interest = known.RNA_modification)

# Extract fractions, Fischer test results, and variants overlapping genes of interest
fraction.RNA_modification <- analysis.RNA_modification[[1]]
test.RNA_modification <- analysis.RNA_modification[[2]]
sig.variants.RNA_mod <- analysis.RNA_modification[[3]]
    
# Add pvalue to the table
all.enrichments['RNA_Modification', 'p'] <- test.RNA_modification$p.value

# Save table of variants overlapping genes of interest
write.table(sig.variants.RNA_mod, file = paste(dir.output.green.SNV, "Variants_overlapping_RNA_modification_genes", ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA", quote = F)

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_RNA_Modification_Genes_Near_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.RNA_modification$Fraction ~ fraction.RNA_modification$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs near RNA modification gene (GO:0009451)',
            ylim = c(0, 0.10),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, .09, paste('p = ', formatC(test.RNA_modification$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.08,
             x1 = 1.88, y1 = 0.08)
    
dev.off()





# CHECK ENRICHMENT OF BRAVO 2024 REGULATORS

# Load files with regulatory genes and combine into 1 df
Bravo_regulators.a <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bravo_2024_GEUVADIS_QTLs/eQTL_Trios_Aggregate.txt', header = TRUE, sep = "\t")
Bravo_regulators.b <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bravo_2024_GEUVADIS_QTLs/eQTL_Trios_Intronic.txt', header = TRUE, sep = "\t")
Bravo_regulators.c <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bravo_2024_GEUVADIS_QTLs/eQTL_Trios_Intergenic_Nearby.txt', header = TRUE, sep = "\t")
Bravo_regulators.d <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bravo_2024_GEUVADIS_QTLs/eQTL_Trios_Intergenic_Distal.txt', header = TRUE, sep = "\t")
Bravo_regulators.e <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bravo_2024_GEUVADIS_QTLs/eQTL_Trios_Exon.txt', header = TRUE, sep = "\t")
Bravo_regulators <- rbind(Bravo_regulators.a, Bravo_regulators.b, Bravo_regulators.c, Bravo_regulators.d, Bravo_regulators.e)

    # Find indices for entries without a gene symbol
    indices.no.symbol <- which(Bravo_regulators$symbol == "")
    
    # Add Ensembl gene to gene symbol column
    Bravo_regulators[indices.no.symbol, 'symbol'] <- Bravo_regulators[indices.no.symbol, 'gene']

    # Generate the gene list
    Bravo_regulators <- Bravo_regulators$symbol
    
    # Only keep unique entries
    Bravo_regulators <- unique(Bravo_regulators)

# Calculate fractions with annotation and Fischer test results
analysis.Bravo_regulators <- Gene_overlap_enrichment_sig_vs_background(annotated_background_variants = all.snvs[, c('RsID', 'Nearby_Genes')], 
                                                                       annotated_significant_variants = sig.snvs.greenlist[, c('SNP', 'Nearby_Genes')], 
                                                                       genes_of_interest = Bravo_regulators)

# Extract fractions, Fischer test results, and variants overlapping genes of interest
fraction.Bravo_regulators <- analysis.Bravo_regulators[[1]]
test.Bravo_regulators <- analysis.Bravo_regulators[[2]]
sig.variants.Bravo_reg <- analysis.Bravo_regulators[[3]]
    
# Add pvalue to the table
all.enrichments['Bravo_Regulators', 'p'] <- test.Bravo_regulators$p.value

# Save table of variants overlapping genes of interest
write.table(sig.variants.Bravo_reg, file = paste(dir.output.green.SNV, "Variants_overlapping_Bravo_Regulators", ".txt", sep=""), row.names = FALSE, col.names = TRUE, sep = "\t", na = "NA", quote = F)

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_Bravo_Regulators_Near_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.Bravo_regulators$Fraction ~ fraction.Bravo_regulators$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs near Bravo 2024 Regulators',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, .19, paste('p = ', formatC(test.Bravo_regulators$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF ALUJ REPEATS

# Modify significant SNP list (blacklist filtered) to have 1 TE per SNV per row
sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(sig.snvs.greenlist, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Modify background SNP list to have 1 gene per SNV per row
background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(all.snvs, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
  
# Calculate proportions in the significant SNP list
    
    # Find index overlap between sig snvs and repeats
    sig.modified1.overlap <- grep('AluJ', sig.snvs.modified1$Repeat, fixed = FALSE)

    # Extract overlap
    sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
    
    # collect repeats for shared snps, if there are any
    sig.modified1.overlap <- aggregate(Repeat ~ SNP, sig.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    sig.in.repeats <- nrow(sig.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # Find index overlap between background snvs and repeats
    background.modified1.overlap <- grep('AluJ', background.snvs.modified1$Repeat, fixed = FALSE)
    
    # Extract overlap
    background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
    
    # MANUALLY CHECK names of unique, overlapping repeats
    unique.repeat.names <- unique(background.modified1.overlap$Repeat)
    
    # collect repeats for shared snps, if there are any
    background.modified1.overlap <- aggregate(Repeat ~ RsID, background.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    background.in.repeats <- nrow(background.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['AluJ', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_AluJ_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping AluJ Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF ALUS REPEATS

# Modify significant SNP list (blacklist filtered) to have 1 TE per SNV per row
sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(sig.snvs.greenlist, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Modify background SNP list to have 1 gene per SNV per row
background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(all.snvs, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
  
# Calculate proportions in the significant SNP list
    
    # Find index overlap between sig snvs and repeats
    sig.modified1.overlap <- grep('AluS', sig.snvs.modified1$Repeat, fixed = FALSE)

    # Extract overlap
    sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
    
    # collect repeats for shared snps, if there are any
    sig.modified1.overlap <- aggregate(Repeat ~ SNP, sig.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    sig.in.repeats <- nrow(sig.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # Find index overlap between background snvs and repeats
    background.modified1.overlap <- grep('AluS', background.snvs.modified1$Repeat, fixed = FALSE)
    
    # Extract overlap
    background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
    
    # MANUALLY CHECK names of unique, overlapping repeats
    unique.repeat.names <- unique(background.modified1.overlap$Repeat)
    
    # collect repeats for shared snps, if there are any
    background.modified1.overlap <- aggregate(Repeat ~ RsID, background.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    background.in.repeats <- nrow(background.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['AluS', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_AluS_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping AluS Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF ALUY REPEATS

# Modify significant SNP list (blacklist filtered) to have 1 TE per SNV per row
sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(sig.snvs.greenlist, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Modify background SNP list to have 1 gene per SNV per row
background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(all.snvs, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
  
# Calculate proportions in the significant SNP list
    
    # Find index overlap between sig snvs and repeats
    sig.modified1.overlap <- grep('AluY', sig.snvs.modified1$Repeat, fixed = FALSE)

    # Extract overlap
    sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
    
    # collect repeats for shared snps, if there are any
    sig.modified1.overlap <- aggregate(Repeat ~ SNP, sig.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    sig.in.repeats <- nrow(sig.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # Find index overlap between background snvs and repeats
    background.modified1.overlap <- grep('AluY', background.snvs.modified1$Repeat, fixed = FALSE)
    
    # Extract overlap
    background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
    
    # MANUALLY CHECK names of unique, overlapping repeats
    unique.repeat.names <- unique(background.modified1.overlap$Repeat)
    
    # collect repeats for shared snps, if there are any
    background.modified1.overlap <- aggregate(Repeat ~ RsID, background.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    background.in.repeats <- nrow(background.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['AluY', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_AluY_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping AluY Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.183,
             x1 = 1.88, y1 = 0.183)
    
dev.off()





# CHECK ENRICHMENT OF L1M REPEATS

# Modify significant SNP list (blacklist filtered) to have 1 TE per SNV per row
sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(sig.snvs.greenlist, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Modify background SNP list to have 1 gene per SNV per row
background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(all.snvs, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
  
# Calculate proportions in the significant SNP list
    
    # Find index overlap between sig snvs and repeats
    sig.modified1.overlap <- grepl('L1M', sig.snvs.modified1$Repeat, fixed = FALSE) & !grepl('HAL1', sig.snvs.modified1$Repeat)

    # Extract overlap
    sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
    
    # collect repeats for shared snps, if there are any
    sig.modified1.overlap <- aggregate(Repeat ~ SNP, sig.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    sig.in.repeats <- nrow(sig.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # Find index overlap between background snvs and repeats
    background.modified1.overlap <- grepl('L1M', background.snvs.modified1$Repeat, fixed = FALSE) & !grepl('HAL1', background.snvs.modified1$Repeat)
    
    # Extract overlap
    background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
    
    # MANUALLY CHECK names of unique, overlapping repeats
    unique.repeat.names <- unique(background.modified1.overlap$Repeat)
    
    # collect repeats for shared snps, if there are any
    background.modified1.overlap <- aggregate(Repeat ~ RsID, background.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    background.in.repeats <- nrow(background.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['L1M', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_L1M_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping L1M Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()












# CHECK ENRICHMENT OF L1P REPEATS

# Modify significant SNP list (blacklist filtered) to have 1 TE per SNV per row
sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(sig.snvs.greenlist, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Modify background SNP list to have 1 gene per SNV per row
background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(all.snvs, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
  
# Calculate proportions in the significant SNP list
    
    # Find index overlap between sig snvs and repeats
    sig.modified1.overlap <- grepl('L1P', sig.snvs.modified1$Repeat, fixed = FALSE) & !grepl('L1PA', sig.snvs.modified1$Repeat)

    # Extract overlap
    sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
    
    # collect repeats for shared snps, if there are any
    sig.modified1.overlap <- aggregate(Repeat ~ SNP, sig.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    sig.in.repeats <- nrow(sig.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # Find index overlap between background snvs and repeats
    background.modified1.overlap <- grepl('L1P', background.snvs.modified1$Repeat, fixed = FALSE) & !grepl('L1PA', background.snvs.modified1$Repeat)
    
    # Extract overlap
    background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
    
    # MANUALLY CHECK names of unique, overlapping repeats
    unique.repeat.names <- unique(background.modified1.overlap$Repeat)
    
    # collect repeats for shared snps, if there are any
    background.modified1.overlap <- aggregate(Repeat ~ RsID, background.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    background.in.repeats <- nrow(background.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['L1P', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_L1P_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping L1P Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF L1PA REPEATS

# Modify significant SNP list (blacklist filtered) to have 1 TE per SNV per row
sig.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(sig.snvs.greenlist, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Modify background SNP list to have 1 gene per SNV per row
background.snvs.modified1 <- as.data.frame(splitstackshape::cSplit(all.snvs, splitCols = "Repeat", sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
  
# Calculate proportions in the significant SNP list
    
    # Find index overlap between sig snvs and repeats
    sig.modified1.overlap <- grep('L1PA|L1HS', sig.snvs.modified1$Repeat, fixed = FALSE)

    # Extract overlap
    sig.modified1.overlap <- sig.snvs.modified1[sig.modified1.overlap, ]
    
    # collect repeats for shared snps, if there are any
    sig.modified1.overlap <- aggregate(Repeat ~ SNP, sig.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    sig.in.repeats <- nrow(sig.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # Find index overlap between background snvs and repeats
    background.modified1.overlap <- grep('L1PA|L1HS', background.snvs.modified1$Repeat, fixed = FALSE)
    
    # Extract overlap
    background.modified1.overlap <- background.snvs.modified1[background.modified1.overlap, ]
    
    # MANUALLY CHECK names of unique, overlapping repeats
    unique.repeat.names <- unique(background.modified1.overlap$Repeat)
    
    # collect repeats for shared snps, if there are any
    background.modified1.overlap <- aggregate(Repeat ~ RsID, background.modified1.overlap, FUN = toString)
    
    # number of unique snps overlapping repeats
    background.in.repeats <- nrow(background.modified1.overlap)
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['L1PA', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_L1PA_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping L1PA/L1HS Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF Non-Intact Full Length L1s (L1Base2)
  
# Calculate proportions in the significant SNP list
    
    # number of unique snps overlapping non-intact fl L1
    sig.in.repeats <- length(which(sig.snvs.greenlist$Nonintact_Full_L1 == 'Yes'))
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # number of unique snps overlapping repeats
    background.in.repeats <- length(which(all.snvs$Nonintact_Full_L1 == 'Yes'))
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['Nonintact-Full_L1', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_L1Base2_Nonintact_Full_L1_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping Nonintact, Full-Length L1 Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF Full Length L1s (L1Base2)
  
# Calculate proportions in the significant SNP list
    
    # number of unique snps overlapping non-intact fl L1
    sig.in.repeats <- length(which(sig.snvs.greenlist$Intact_Full_L1 == 'Yes'))
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # number of unique snps overlapping repeats
    background.in.repeats <- length(which(all.snvs$Intact_Full_L1 == 'Yes'))
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['Intact-Full_L1', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_L1Base2_Intact_Full_L1_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping Intact, Full-Length L1 Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF ORF2-Intact L1s (L1Base2)
  
# Calculate proportions in the significant SNP list
    
    # number of unique snps overlapping non-intact fl L1
    sig.in.repeats <- length(which(sig.snvs.greenlist$Intact_ORF2 == 'Yes'))
    
    # number of unique snps not overlapping repeats
    sig.notin.repeats <- nrow(sig.snvs.greenlist) - sig.in.repeats
    
    # fraction of snps overlapping repeats
    sig.fraction.repeats <- sig.in.repeats / (sig.in.repeats + sig.notin.repeats)

# Calculate proportions in the background SNP list
    
    # number of unique snps overlapping repeats
    background.in.repeats <- length(which(all.snvs$Intact_ORF2 == 'Yes'))
    
    # number of unique snps not overlapping repeats
    background.notin.repeats <- nrow(all.snvs) - background.in.repeats
    
    # fraction of snps overlapping repeats
    background.fraction.repeats <- background.in.repeats / (background.in.repeats + background.notin.repeats)
    
# Collect fractions in a df to use for plot
fraction.repeats <- data.frame(Fraction = c(sig.fraction.repeats, background.fraction.repeats), Group = c('Significant', 'Background'))
    
# Generate a contingency table as a matrix to run stats
my.contingency.repeats <- as.matrix(data.frame(Not_Repeat = c(sig.notin.repeats, background.notin.repeats), Repeat = c(sig.in.repeats, background.in.repeats)))

# Run Fischer's Exact Test
test.repeats <- fisher.test(my.contingency.repeats)
    
# Add pvalue to the table
all.enrichments['Intact_ORF2_L1', 'p'] <- test.repeats$p.value

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_L1Base2_Intact_ORF2_L1_overlapping_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.repeats$Fraction ~ fraction.repeats$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs overlapping Intact ORF2 L1 Repeats',
            ylim = c(0, 0.20),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.19, paste('p = ', formatC(test.repeats$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.18,
             x1 = 1.88, y1 = 0.18)
    
dev.off()





# CHECK ENRICHMENT OF ENCODE cCREs

# Calculate fractions with annotation and Fischer test results
analysis.cCREs <- Annotation_enrichment_sig_vs_background(annotated_background_variants = all.snvs, 
                                                                   annotated_significant_variants = sig.snvs.greenlist, 
                                                                   annotation_column_name = 'ENCODE_cCRE', 
                                                                   annotation_label = 'Yes')
    
# Extract fractions and Fischer test results
fraction.cCRE <- analysis.cCREs[[1]]
test.cCRE <- analysis.cCREs[[2]]
  
# Add pvalue to the table
all.enrichments['ENCODE_cCREs', 'p'] <- test.cCRE$p.value  

# Generate bar plots comparing frequency of blacklisted SNVs
pdf(paste(dir.output.green.SNV, "Barplot_Frequency_ENCODE_cCREs_Greenlist_SNVs", ".pdf", sep=""), width = 6, height = 6)
par(mar=c(6,4,4,1)+.1)

    # Plot proportion mediated as a barplot
    barplot(fraction.cCRE$Fraction ~ fraction.cCRE$Group,
            ylab = "Fraction of SNVs",
            main = 'SNVs in ENCODE Registry v4 cCRE',
            ylim = c(0, 0.30),
            xlab = "",
            las = 2,
            cex.names = 1
            )
    
    # add p-values
    text(1.25, 0.28, paste('p = ', formatC(test.cCRE$p.value, format = 'E', digits = 2), sep = ''), cex = 1)
    
    # Draw line between the two grouos
    segments(x0 = .7, y0 = 0.27,
             x1 = 1.88, y1 = 0.27)
    
dev.off()





# Calculate final FDR values
all.enrichments[, 'FDR'] <- p.adjust(all.enrichments[, 'p'], method = 'fdr')

    # Order by pvalue
    all.enrichments <- all.enrichments[order(all.enrichments$p, decreasing = FALSE), ]

# Save results
write.table(all.enrichments, file = paste(dir.output.green.SNV, "Greenlist_SNV_Enrichment_pvalues_FDR", ".txt", sep=""), row.names = TRUE, col.names = NA, sep = "\t", na = "NA", quote = F)



    

# Clean the environment
rm(list=ls())       
        

