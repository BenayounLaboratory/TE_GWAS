# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(data.table) # for fread and fsave

# Specify directory for annotated SNVs
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/'
dir.output.EXTRA.BED <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/BED_FILES_WITH_WINDOW/'


# CLEANUP GREAT ANNOTATIONS ------------------------------------------------------------------


# Specify the directory for GREAT annotation files
dir.GREAT <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/'

# Collect individual file paths in a list.
GREAT_list <- list.files(dir.GREAT, "\\GREAT_Set", recursive=FALSE, full.names=TRUE)
  
# Create new dataframe to hold all annotations
GREAT_annotations <- data.frame()

# Start a counter
GREAT_counter <- 0

# Loop over each parquet file path
for(GREAT_file in GREAT_list) {
  
    # Update the GREAT file counter
    GREAT_counter <- GREAT_counter + 1
    
    # Read in the ith GREAT file
    ith_GREAT_data <- fread(file = GREAT_file, header = FALSE, sep = '\t', data.table = FALSE)
    
    # Define or extend a df with ith data
    if (GREAT_counter == 1) {
      
      # If this is the first GREAT file, define a new df to hold the data. 
      GREAT_annotations <- ith_GREAT_data
      
    } else {
      
      # If this is NOT the first GREAT file, add the data to the existing df.
      GREAT_annotations <- rbind(GREAT_annotations, ith_GREAT_data)
      
    }
  
} # END FOR LOOP

# Remove temporary variables
rm(ith_GREAT_data)

# Update colnames
colnames(GREAT_annotations) <- c('RsID', 'Nearby_Genes')

# Check if SNPs are unique
nrow(GREAT_annotations) == sum(!duplicated(GREAT_annotations$RsID))

# Remove gene distances
GREAT_annotations$Nearby_Genes <- gsub("\\s*\\([^\\)]+\\)","", as.character(GREAT_annotations$Nearby_Genes))

# Save the csv
fwrite(GREAT_annotations, paste(dir.GREAT, "GREAT_Annotations_Cleaned_Up.csv", sep = ''))



# APPLY ANNOTATIONS TO SNVS ------------------------------------------------------------------
        
        
# PREPARE ANNOTATIONS

# Load SNPEff results (FOR SIGNIFICANT variants only, not for background variants) (I'm only annotating phastcon scores)
annot.SNPEff <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNPEff/READABLE_ANNOTATIONS.csv', header = TRUE, sep = '\t')

    # Remove duplicate RsIDs
    annot.SNPEff <- annot.SNPEff[!duplicated(annot.SNPEff$ID), ]


# Load ENCODE Blacklist annotation files
annot.Blacklist.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_ENCODE_Blacklist.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.Blacklist.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_ENCODE_Blacklist.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.Blacklist <- rbind(annot.Blacklist.SNV, annot.Blacklist.SV)

    # Only keep necessary columns
    annot.Blacklist <- annot.Blacklist[, c('V4', 'V8')]
    
    # Update colnames
    colnames(annot.Blacklist) <- c('RsID', 'Blacklist_Label')
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(annot.Blacklist$RsID)) # 0
    
    # Assign RsID to rownames
    rownames(annot.Blacklist) <- annot.Blacklist$RsID
    
    # remove unneeded variables
    rm(annot.Blacklist.SNV, annot.Blacklist.SV)
    
    
# Load SV Hotspot (Lin 2019) annotation files
annot.SVHotspot.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_SV-Hotspots_Lin-2019.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SVHotspot.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_SV-Hotspots_Lin-2019.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SVHotspot <- rbind(annot.SVHotspot.SNV, annot.SVHotspot.SV)

    # Only keep necessary columns
    annot.SVHotspot <- annot.SVHotspot[, c('V4', 'V8')]
    
    # Update colnames
    colnames(annot.SVHotspot) <- c('RsID', 'SV_Hotspot_Label')
    
    # Assign a different value for the the SV_hotspot_label column
    annot.SVHotspot$SV_Hotspot_Label <- 'Lin_2019'
      
    # Check whether there are duplicate RsIDs (i.e. overlap more than 1 region)
    sum(duplicated(annot.SVHotspot$RsID)) # 2
    
    # Remove duplicate RsIDs
    annot.SVHotspot <- annot.SVHotspot[!duplicated(annot.SVHotspot$RsID), ]
    
    # Assign RsID to rownames
    rownames(annot.SVHotspot) <- annot.SVHotspot$RsID
    
    # remove unneeded variables
    rm(annot.SVHotspot.SNV, annot.SVHotspot.SV)
    
    
# Load Segmental Duplications (Bailey 2001) annotation files
annot.SegmentalDup.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_Segmental_Duplications_Bailey-2001.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SegmentalDup.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_Segmental_Duplications_Bailey-2001.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SegmentalDup <- rbind(annot.SegmentalDup.SNV, annot.SegmentalDup.SV)

    # Only keep necessary columns
    annot.SegmentalDup <- annot.SegmentalDup[, c('V4', 'V8')]
    
    # Update colnames
    colnames(annot.SegmentalDup) <- c('RsID', 'Segmental_Duplication')
    
    # Assign a different value for the the label column
    annot.SegmentalDup$Segmental_Duplication <- 'Bailey_2001'
      
    # Check whether there are duplicate RsIDs (i.e. overlap more than 1 region)
    sum(duplicated(annot.SegmentalDup$RsID)) # 623240
    
    # Remove duplicate RsIDs
    annot.SegmentalDup <- annot.SegmentalDup[!duplicated(annot.SegmentalDup$RsID), ]
    
    # Assign RsID to rownames
    rownames(annot.SegmentalDup) <- annot.SegmentalDup$RsID
    
    # remove unneeded variables
    rm(annot.SegmentalDup.SNV, annot.SegmentalDup.SV)
            
    
# Load GREAT annotation file (SNVs and SVs)
GREAT_annotations <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/GREAT_Annotations_Cleaned_Up.csv', header = TRUE, sep = ',', data.table = FALSE)
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(GREAT_annotations$RsID)) # 0
    
    # Assign RsID to rownames
    rownames(GREAT_annotations) <- GREAT_annotations$RsID
            
    
# Load Bravo 2023 GEUVADIS LCL cis-eQTLs (only for SNVs)
Bravo_eQTLs <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bravo_2024_GEUVADIS_QTLs/Gene_cis-eQTLs.txt', header = TRUE, sep = '\t', data.table = FALSE)
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(Bravo_eQTLs$RsID)) # 297599
    
    # Only keep necessary columns
    Bravo_eQTLs <- Bravo_eQTLs[, c('RsID', 'gene')]
    
    # Update colnames
    colnames(Bravo_eQTLs) <- c('RsID', 'cis_eQTL')
      
    # Aggregate cis_eQTLs by RsID
    Bravo_eQTLs <- aggregate(cis_eQTL ~ RsID, Bravo_eQTLs, FUN = toString)
    
        # Update the rownames
        rownames(Bravo_eQTLs) <- Bravo_eQTLs$RsID
        
    # Check whether there are duplicate RsIDs
    sum(duplicated(Bravo_eQTLs$RsID)) # 0
    
    
# Load REPEATMASKER annotation files
annot.Repeatmasker.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_Repeatmasker.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.Repeatmasker.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_Repeatmasker.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.Repeatmasker <- rbind(annot.Repeatmasker.SNV, annot.Repeatmasker.SV)

    # Add column with TE lengths
    annot.Repeatmasker$Repeat_length <- abs(annot.Repeatmasker$V7 - annot.Repeatmasker$V6)
    
    # Only keep necessary columns
    annot.Repeatmasker <- annot.Repeatmasker[, c('V4', 'V8', 'Repeat_length')]
    
    # Update colnames
    colnames(annot.Repeatmasker) <- c('RsID', 'Repeat', 'Repeat_length')
      
    # Aggregate TE names by RsID
    aggregate.TEs <- aggregate(Repeat ~ RsID, annot.Repeatmasker, FUN = toString)
    
        # Update the rownames
        rownames(aggregate.TEs) <- aggregate.TEs$RsID
        
    # Aggregate TE lengths by RsID
    aggregate.TE_lengths <- aggregate(Repeat_length ~ RsID, annot.Repeatmasker, FUN = toString)
    
        # Update the rownames
        rownames(aggregate.TE_lengths) <- aggregate.TE_lengths$RsID
    
    # Replace repeatmasker annotations with aggregate results
    annot.Repeatmasker <- aggregate.TEs
    annot.Repeatmasker$Repeat_length <- aggregate.TE_lengths$Repeat_length
    
    # Check whether there are duplicate RsIDs
    sum(duplicated(annot.Repeatmasker$RsID)) # 0
   
    # remove unneeded variables
    rm(annot.Repeatmasker.SNV, annot.Repeatmasker.SV, aggregate.TEs, aggregate.TE_lengths)
        

# Load L1Base2 full length, intact L1 overlaps
annot.L1Base.full.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_L1Base2_Full_Length_L1s.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.L1Base.full.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_L1Base2_Full_Length_L1s.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.L1Base.full <- rbind(annot.L1Base.full.SNV, annot.L1Base.full.SV)

    # Add a label
    annot.L1Base.full$Intact_Full_L1 <- 'Yes'
      
    # Only keep necessary columns
    annot.L1Base.full <- annot.L1Base.full[, c('V4', 'Intact_Full_L1')]
    
    # Update colnames
    colnames(annot.L1Base.full) <- c('RsID', 'Intact_Full_L1')
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(annot.L1Base.full$RsID)) # 0
    
    # Assign RsID to rownames
    rownames(annot.L1Base.full) <- annot.L1Base.full$RsID    
    
    
# Load L1Base2 ORF2-intact L1 overlaps
annot.L1Base.ORF2.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_L1Base2_Intact_ORF2_L1.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.L1Base.ORF2.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_L1Base2_Intact_ORF2_L1.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.L1Base.ORF2 <- rbind(annot.L1Base.ORF2.SNV, annot.L1Base.ORF2.SV)

    # Add a label
    annot.L1Base.ORF2$Intact_ORF2 <- 'Yes'
      
    # Only keep necessary columns
    annot.L1Base.ORF2 <- annot.L1Base.ORF2[, c('V4', 'Intact_ORF2')]
    
    # Update colnames
    colnames(annot.L1Base.ORF2) <- c('RsID', 'Intact_ORF2')
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(annot.L1Base.ORF2$RsID)) # 0
    
    # Assign RsID to rownames
    rownames(annot.L1Base.ORF2) <- annot.L1Base.ORF2$RsID    
    
    
# Load L1Base2 full length, nonintact L1 overlaps
annot.L1Base.nonintact.SNV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNV_Variants_Intersecting_L1Base2_NonIntact_Full_Length_L1s.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.L1Base.nonintact.SV <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SV_Variants_Intersecting_L1Base2_NonIntact_Full_Length_L1s.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.L1Base.nonintact <- rbind(annot.L1Base.nonintact.SNV, annot.L1Base.nonintact.SV)

    # Add a label
    annot.L1Base.nonintact$Nonintact_Full_L1 <- 'Yes'
      
    # Only keep necessary columns
    annot.L1Base.nonintact <- annot.L1Base.nonintact[, c('V4', 'Nonintact_Full_L1')]
    
    # Update colnames
    colnames(annot.L1Base.nonintact) <- c('RsID', 'Nonintact_Full_L1')
    
    # Remove duplicated SNVs
    annot.L1Base.nonintact <- annot.L1Base.nonintact[!duplicated(annot.L1Base.nonintact), ]
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(annot.L1Base.nonintact$RsID)) # 0
    
    # Assign RsID to rownames
    rownames(annot.L1Base.nonintact) <- annot.L1Base.nonintact$RsID    
    
    
# Load ENCODE cCRE Annotations
annot.SNV.cCRE <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SNV_Variants_Intersecting_ENCODE_cCREs.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SV.cCRE <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SV_Variants_Intersecting_ENCODE_cCREs.txt', header = FALSE, sep = '\t', data.table = FALSE)

annot.SNV.promoter <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SNV_Variants_Intersecting_ENCODE_Candidate_Promoters.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SV.promoter <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SV_Variants_Intersecting_ENCODE_Candidate_Promoters.txt', header = FALSE, sep = '\t', data.table = FALSE)

annot.SNV.enhancer <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SNV_Variants_Intersecting_ENCODE_Candidate_Enhancers.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SV.enhancer <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SV_Variants_Intersecting_ENCODE_Candidate_Enhancers.txt', header = FALSE, sep = '\t', data.table = FALSE)

annot.SNV.CTCF <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SNV_Variants_Intersecting_ENCODE_CTCF_Bound.txt', header = FALSE, sep = '\t', data.table = FALSE)
annot.SV.CTCF <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations/SV_Variants_Intersecting_ENCODE_CTCF_Bound.txt', header = FALSE, sep = '\t', data.table = FALSE)

    # Combine SNVs and SVs
    annot.cCRE <- rbind(annot.SNV.cCRE, annot.SV.cCRE)
    annot.promoter <- rbind(annot.SNV.promoter, annot.SV.promoter)
    annot.enhancer <- rbind(annot.SNV.enhancer, annot.SV.enhancer)
    annot.CTCF <- rbind(annot.SNV.CTCF, annot.SV.CTCF)

    # Add a label
    annot.cCRE$ENCODE_cCRE <- 'Yes'
    annot.promoter$ENCODE_cPromoter <- 'Yes'
    annot.enhancer$ENCODE_cEnhancer <- 'Yes'
    annot.CTCF$ENCODE_CTCF <- 'Yes'
      
    # Only keep necessary columns
    annot.cCRE <- annot.cCRE[, c('V4', 'ENCODE_cCRE')]
    annot.promoter <- annot.promoter[, c('V4', 'ENCODE_cPromoter')]
    annot.enhancer <- annot.enhancer[, c('V4', 'ENCODE_cEnhancer')]
    annot.CTCF <- annot.CTCF[, c('V4', 'ENCODE_CTCF')]
    
    # Update colnames
    colnames(annot.cCRE)[1] <- c('RsID')
    colnames(annot.promoter)[1] <- c('RsID')
    colnames(annot.enhancer)[1] <- c('RsID')
    colnames(annot.CTCF)[1] <- c('RsID')
    
    # Remove duplicated variants
    annot.cCRE <- annot.cCRE[!duplicated(annot.cCRE), ]
    annot.promoter <- annot.promoter[!duplicated(annot.promoter), ]
    annot.enhancer <- annot.enhancer[!duplicated(annot.enhancer), ]
    annot.CTCF <- annot.CTCF[!duplicated(annot.CTCF), ]
      
    # Check whether there are duplicate RsIDs
    sum(duplicated(annot.cCRE$RsID)) # 0
    sum(duplicated(annot.promoter$RsID)) # 0
    sum(duplicated(annot.enhancer$RsID)) # 0
    sum(duplicated(annot.CTCF$RsID)) # 0
    
    # Assign RsID to rownames
    rownames(annot.cCRE) <- annot.cCRE$RsID 
    rownames(annot.promoter) <- annot.promoter$RsID 
    rownames(annot.enhancer) <- annot.enhancer$RsID 
    rownames(annot.CTCF) <- annot.CTCF$RsID 
    
      
    
    
    
    
    
    
    
    
# APPLY ANNOTATIONS TO ALL GWAS SNPS/SVs
 
# Load list of SNVs only
list.SNVs <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_IDs_SNVs_only.txt', header = FALSE, sep = "\t")
   
# Load list of SVs only
list.SVs <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_IDs_SVs_only.txt', header = FALSE, sep = "\t")

# Load the file to annotate
all.variants <- fread(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_Variant_Map_AFR_SNVs_and_SVs.txt', header = TRUE, sep = "\t", data.table = FALSE)
   
    # Remove unnecessary columns
    all.variants <- all.variants[, -c(1)]
    
    # Check whether there are duplicate RsIDs
    sum(duplicated(all.variants$RsID)) # 0
    
    # Update rownames with RsIDs
    rownames(all.variants) <- all.variants$RsID
        
# Add columns to the variant table to hold the annotations
all.variants$Blacklist_Label <- NA
all.variants$SV_Hotspot_Label <- NA
all.variants$Segmental_Duplication <- NA
all.variants$Nearby_Genes <- NA
all.variants$Bravo_eQTLs <- NA
all.variants$Repeat <- NA 
all.variants$Repeat_length <- NA 
all.variants$Intact_Full_L1 <- NA 
all.variants$Intact_ORF2 <- NA 
all.variants$Nonintact_Full_L1 <- NA 
all.variants$ENCODE_cCRE <- NA 
all.variants$ENCODE_cPromoter <- NA 
all.variants$ENCODE_cEnhancer <- NA 
all.variants$ENCODE_CTCF <- NA 

# Find indices for variants that overlap with each annotation
SNPs.Blacklist <-    which(all.variants$RsID %in% annot.Blacklist$RsID)
SNPs.SV_Hotspots <-    which(all.variants$RsID %in% annot.SVHotspot$RsID)
SNPs.SegmentalDup <-    which(all.variants$RsID %in% annot.SegmentalDup$RsID)
SNPs.GREAT <-        which(all.variants$RsID %in% GREAT_annotations$RsID)
SNPs.Bravo_eQTLs <-        which(all.variants$RsID %in% Bravo_eQTLs$RsID)
SNPs.Repeatmasker <- which(all.variants$RsID %in% annot.Repeatmasker$RsID)        
SNPs.L1.full <- which(all.variants$RsID %in% annot.L1Base.full$RsID)         
SNPs.L1.ORF2 <- which(all.variants$RsID %in% annot.L1Base.ORF2$RsID)         
SNPs.nonintact <- which(all.variants$RsID %in% annot.L1Base.nonintact$RsID)      
SNPs.cCRE <- which(all.variants$RsID %in% annot.cCRE$RsID)    
SNPs.promoter <- which(all.variants$RsID %in% annot.promoter$RsID)    
SNPs.enhancer <- which(all.variants$RsID %in% annot.enhancer$RsID)    
SNPs.CTCF <- which(all.variants$RsID %in% annot.CTCF$RsID)    
            
# Define the RsIDs corresponding to indices
SNPs.Blacklist <-     all.variants[SNPs.Blacklist, 'RsID']        
SNPs.SV_Hotspots <-     all.variants[SNPs.SV_Hotspots, 'RsID']        
SNPs.SegmentalDup <-     all.variants[SNPs.SegmentalDup, 'RsID']        
SNPs.GREAT <-         all.variants[SNPs.GREAT, 'RsID']   
SNPs.Bravo_eQTLs <-         all.variants[SNPs.Bravo_eQTLs, 'RsID']   
SNPs.Repeatmasker <-  all.variants[SNPs.Repeatmasker, 'RsID']        
SNPs.L1.full <- all.variants[SNPs.L1.full, 'RsID']         
SNPs.L1.ORF2 <- all.variants[SNPs.L1.ORF2, 'RsID']            
SNPs.nonintact <- all.variants[SNPs.nonintact, 'RsID']         
SNPs.cCRE <- all.variants[SNPs.cCRE, 'RsID']  
SNPs.promoter <- all.variants[SNPs.promoter, 'RsID']  
SNPs.enhancer <- all.variants[SNPs.enhancer, 'RsID']  
SNPs.CTCF <- all.variants[SNPs.CTCF, 'RsID']  

# Update the variant table
all.variants[SNPs.Blacklist, 'Blacklist_Label'] <- annot.Blacklist[SNPs.Blacklist, 'Blacklist_Label']
all.variants[SNPs.SV_Hotspots, 'SV_Hotspot_Label'] <- annot.SVHotspot[SNPs.SV_Hotspots, 'SV_Hotspot_Label']
all.variants[SNPs.SegmentalDup, 'Segmental_Duplication'] <- annot.SegmentalDup[SNPs.SegmentalDup, 'Segmental_Duplication']
all.variants[SNPs.GREAT, 'Nearby_Genes'] <- GREAT_annotations[SNPs.GREAT, 'Nearby_Genes']
all.variants[SNPs.Bravo_eQTLs, 'Bravo_eQTLs'] <- Bravo_eQTLs[SNPs.Bravo_eQTLs, 'cis_eQTL']
all.variants[SNPs.Repeatmasker, c('Repeat', 'Repeat_length')] <- annot.Repeatmasker[SNPs.Repeatmasker, c('Repeat', 'Repeat_length')]
all.variants[SNPs.L1.full, c('Intact_Full_L1')] <- annot.L1Base.full[SNPs.L1.full, c('Intact_Full_L1')]
all.variants[SNPs.L1.ORF2, c('Intact_ORF2')] <- annot.L1Base.ORF2[SNPs.L1.ORF2, c('Intact_ORF2')]
all.variants[SNPs.nonintact, c('Nonintact_Full_L1')] <- annot.L1Base.nonintact[SNPs.nonintact, c('Nonintact_Full_L1')]
all.variants[SNPs.cCRE, c('ENCODE_cCRE')] <- annot.cCRE[SNPs.cCRE, c('ENCODE_cCRE')]
all.variants[SNPs.promoter, c('ENCODE_cPromoter')] <- annot.promoter[SNPs.promoter, c('ENCODE_cPromoter')]
all.variants[SNPs.enhancer, c('ENCODE_cEnhancer')] <- annot.enhancer[SNPs.enhancer, c('ENCODE_cEnhancer')]
all.variants[SNPs.CTCF, c('ENCODE_CTCF')] <- annot.CTCF[SNPs.CTCF, c('ENCODE_CTCF')]

# Separately save annotated SNVs and SVs
fwrite(all.variants[which(all.variants$RsID %in% list.SNVs$V1), ], file = paste(dir.output, "All_Annotated_GWAS_Background_SNVs_only.csv", sep=""), na = 'NA')     
fwrite(all.variants[which(all.variants$RsID %in% list.SVs$V1), ], file = paste(dir.output, "All_Annotated_GWAS_Background_SVs_only.csv", sep=""), na = 'NA')     


    
    
    
    
    
    
    
    
# MAKE GENE LISTS FOR BACKGROUND
        
# Subset SNV and SV variants
genes.background.SNV <- all.variants[which(all.variants$RsID %in% list.SNVs$V1), ]
genes.background.SV <- all.variants[which(all.variants$RsID %in% list.SVs$V1), ]

# only keep variant IDs and nearby genes
genes.background.SNV <- genes.background.SNV[, c('RsID', 'Nearby_Genes')]
genes.background.SV <- genes.background.SV[, c('RsID', 'Nearby_Genes')]

# Modify significant variant lists to have 1 gene per SNV per row
genes.background.SNV <- as.data.frame(splitstackshape::cSplit(genes.background.SNV, splitCols = 'Nearby_Genes', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
genes.background.SV <- as.data.frame(splitstackshape::cSplit(genes.background.SV, splitCols = 'Nearby_Genes', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Remove duplicate genes
genes.background.SNV <- genes.background.SNV[!duplicated(genes.background.SNV$Nearby_Genes), ]
genes.background.SV <- genes.background.SV[!duplicated(genes.background.SV$Nearby_Genes), ]

# Remove NAs
genes.background.SNV <- genes.background.SNV[!is.na(genes.background.SNV$Nearby_Genes), ]
genes.background.SV <- genes.background.SV[!is.na(genes.background.SV$Nearby_Genes), ]

# Remove ENSEMBL genes, only keep gene symbols
genes.background.SNV <- genes.background.SNV[!grepl('ENSG', genes.background.SNV$Nearby_Genes), ]
genes.background.SV <- genes.background.SV[!grepl('ENSG', genes.background.SV$Nearby_Genes), ]

# Save gene lists
write.table(genes.background.SNV[, 'Nearby_Genes', drop = FALSE], file = paste(dir.output, "NEARBY_GENE_LIST_BACKGROUND_SNVs.txt", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
write.table(genes.background.SV[, 'Nearby_Genes', drop = FALSE], file = paste(dir.output, "NEARBY_GENE_LIST_BACKGROUND_SVs.txt", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)          
    
    
    






# APPLY ANNOTATIONS TO SIGNIFICANT GWAS SNVs/SVs  
        
# Load significant SNVs/SVs        
sig.variants <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/Singleton_Meta-Analysis_Significant_Variants.txt', header = TRUE, sep = "")    

    # Assign RsIDs to rownames
    rownames(sig.variants) <- sig.variants$SNP
    
# Load clumping results        
clumped.variants <- read.csv(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/Singleton_Meta-Analysis_Significant_Variants.clumped', header = TRUE, sep = "")    

    # Assign RsIDs to rownames
    rownames(clumped.variants) <- clumped.variants$SNP
        
# Add columns to the significant SNV/SV table to hold the annotations
sig.variants$PhastCons <- NA
sig.variants$Blacklist_Label <- NA
sig.variants$SV_Hotspot_Label <- NA
sig.variants$Segmental_Duplication <- NA
sig.variants$Nearby_Genes <- NA
sig.variants$Bravo_eQTLs <- NA
sig.variants$Repeat <- NA 
sig.variants$Repeat_length <- NA 
sig.variants$Intact_Full_L1 <- NA 
sig.variants$Intact_ORF2 <- NA 
sig.variants$Nonintact_Full_L1 <- NA 
sig.variants$Index_Variant <- 'No' 
sig.variants$Index_Variant_ID <- '' 
sig.variants$ENCODE_cCRE <- NA
sig.variants$ENCODE_cPromoter <- NA
sig.variants$ENCODE_cEnhancer <- NA
sig.variants$ENCODE_CTCF <- NA
            
# Define variants to subset from the entire annotated variant table
SNPs.to.subset <- sig.variants$SNP

    # Update the significant SNV/SV table
    sig.variants[SNPs.to.subset, c('Blacklist_Label', 'SV_Hotspot_Label', 'Segmental_Duplication', 'Nearby_Genes', 'Bravo_eQTLs', 'Repeat', 'Repeat_length', 'Intact_Full_L1', 'Intact_ORF2', 'Nonintact_Full_L1', 'ENCODE_cCRE', 'ENCODE_cPromoter', 'ENCODE_cEnhancer', 'ENCODE_CTCF')] <- all.variants[SNPs.to.subset, c('Blacklist_Label', 'SV_Hotspot_Label', 'Segmental_Duplication', 'Nearby_Genes', 'Bravo_eQTLs', 'Repeat', 'Repeat_length', 'Intact_Full_L1', 'Intact_ORF2', 'Nonintact_Full_L1', 'ENCODE_cCRE', 'ENCODE_cPromoter', 'ENCODE_cEnhancer', 'ENCODE_CTCF')]

    # Add PhastCon scores
    sig.variants[annot.SNPEff$ID, 'PhastCons'] <- annot.SNPEff[, 'PhastCons']
    
# Define variants that are index variants in the clumping analysis
SNPs.index <- which(sig.variants$SNP %in% clumped.variants$SNP)
SNPs.index <- sig.variants[SNPs.index, 'SNP'] 

    # Update the significant SNV/SV table
    sig.variants[SNPs.index, 'Index_Variant'] <- 'Yes' # Temporarily assign yes to all variants

    
# List index SNV each SNP is clumped to, if any

# Loop over each significant SNV
for (ith_SNP in 1:nrow(sig.variants)) {
  
    if (sig.variants[ith_SNP, 'Index_Variant'] == 'Yes') {
      
        # Fill in the cell with the ith SNP
        sig.variants[ith_SNP, 'Index_Variant_ID'] <- sig.variants[ith_SNP, 'SNP']
      
    } else {
      
        # Extract the ith SNP ID
        ith_SNP_ID <- sig.variants[ith_SNP, 'SNP']
        
        # Overlap SNP ID with the clumping results and specify the row for the overlap
        clumping_row <- grep(paste(ith_SNP_ID, '(1)', sep = ""), clumped.variants$SP2, fixed = TRUE)
        
        # Assign the index SNV from the clumping results to the annotation table
        sig.variants[ith_SNP, 'Index_Variant_ID'] <- clumped.variants[clumping_row, 'SNP']
          
    }
  
} # END FOR LOOP over each significant SNV

# Cleanup table, by changing "NONE" annotations in to the nearby gene column to NA
sig.variants[which(sig.variants$Nearby_Genes == "NONE"), 'Nearby_Genes'] <- NA

    # Save ALL results
    fwrite(sig.variants, file = paste(dir.output, "All_Annotated_GWAS_Significant_SNVs_AND_SVs.csv", sep=""), na = 'NA')         
        
# Split significant SNVs/SVs, and also split into blacklist and greenlist tables
    
    # Identify indices for blacklist entries
    Blacklist.indices <- which(!is.na(sig.variants$Blacklist_Label))
    
    # Identify IDs for blacklist entries
    Blacklist.IDs <- sig.variants[Blacklist.indices, 'SNP']
    
    # Identify indices for all variants in LD with blacklisted entries
    Blacklist.LD.indices <- which(sig.variants$Index_Variant_ID %in% Blacklist.IDs)
    
    # Get all blacklist indices, either directly in the blacklist or in LD with the blacklist
    Blacklist.indices.all <- union(Blacklist.indices, Blacklist.LD.indices)
    
    # Split variants into green/blacklist
    variants.blacklist.all <- sig.variants[Blacklist.indices.all, ]
    variants.greenlist.all <- sig.variants[-which(sig.variants$SNP %in% variants.blacklist.all$SNP), ]
    
    # Split blacklists/greenlists into SNVs and SVs
    variants.blacklist.SNV <- variants.blacklist.all[which(variants.blacklist.all$SNP %in% list.SNVs$V1), ]
    variants.blacklist.SV <- variants.blacklist.all[which(variants.blacklist.all$SNP %in% list.SVs$V1), ]
    
    variants.greenlist.SNV <- variants.greenlist.all[which(variants.greenlist.all$SNP %in% list.SNVs$V1), ]
    variants.greenlist.SV <- variants.greenlist.all[which(variants.greenlist.all$SNP %in% list.SVs$V1), ]   
        
    # Save results
    fwrite(variants.blacklist.SNV, file = paste(dir.output, "Annotated_Significant_Variants_Blacklist_SNVs.csv", sep=""))         
    # fwrite(variants.blacklist.SV, file = paste(dir.output, "Annotated_Significant_Variants_Blacklist_SVs.csv", sep="")) # No entries
    fwrite(variants.greenlist.SNV, file = paste(dir.output, "Annotated_Significant_Variants_Greenlist_SNVs.csv", sep=""))     
    fwrite(variants.greenlist.SV, file = paste(dir.output, "Annotated_Significant_Variants_Greenlist_SVs.csv", sep=""))         


    
    
    
    
    
    
    
    
# MAKE GENE LISTS FOR SIGNIFICANT VARIANTS
        
# Modify significant variant lists to have 1 gene per SNV per row
genes.blacklist.SNV <- as.data.frame(splitstackshape::cSplit(variants.blacklist.SNV, splitCols = 'Nearby_Genes', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
genes.greenlist.SNV <- as.data.frame(splitstackshape::cSplit(variants.greenlist.SNV, splitCols = 'Nearby_Genes', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))
genes.greenlist.SV <- as.data.frame(splitstackshape::cSplit(variants.greenlist.SV, splitCols = 'Nearby_Genes', sep = ",", direction = "long", fixed = TRUE, drop = TRUE, stripWhite = TRUE, makeEqual = FALSE))

# Remove duplicate genes
genes.blacklist.SNV <- genes.blacklist.SNV[!duplicated(genes.blacklist.SNV$Nearby_Genes), ]
genes.greenlist.SNV <- genes.greenlist.SNV[!duplicated(genes.greenlist.SNV$Nearby_Genes), ]
genes.greenlist.SV <- genes.greenlist.SV[!duplicated(genes.greenlist.SV$Nearby_Genes), ]

# Remove NAs
genes.blacklist.SNV <- genes.blacklist.SNV[!is.na(genes.blacklist.SNV$Nearby_Genes), ]
genes.greenlist.SNV <- genes.greenlist.SNV[!is.na(genes.greenlist.SNV$Nearby_Genes), ]
genes.greenlist.SV <- genes.greenlist.SV[!is.na(genes.greenlist.SV$Nearby_Genes), ]

# Remove ENSEMBL genes, only keep gene symbols
genes.blacklist.SNV <- genes.blacklist.SNV[!grepl('ENSG', genes.blacklist.SNV$Nearby_Genes), ]
genes.greenlist.SNV <- genes.greenlist.SNV[!grepl('ENSG', genes.greenlist.SNV$Nearby_Genes), ]
genes.greenlist.SV <- genes.greenlist.SV[!grepl('ENSG', genes.greenlist.SV$Nearby_Genes), ]

# Save gene lists
write.table(genes.blacklist.SNV[, 'Nearby_Genes', drop = FALSE], file = paste(dir.output, "NEARBY_GENE_LIST_Significant_Variants_Blacklist_SNVs.txt", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
write.table(genes.greenlist.SNV[, 'Nearby_Genes', drop = FALSE], file = paste(dir.output, "NEARBY_GENE_LIST_Significant_Variants_Greenlist_SNVs.txt", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
write.table(genes.greenlist.SV[, 'Nearby_Genes', drop = FALSE], file = paste(dir.output, "NEARBY_GENE_LIST_Significant_Variants_Greenlist_SVs.txt", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         


    
    
    
    
    
    

    
# MAKE BED FILES 
    
# To make BED files, update the chromosome labels to include 'chr' (needed by pipelines like GREAT)
variants.blacklist.SNV$CHR <- paste('chr', variants.blacklist.SNV$CHR, sep = '')
#variants.blacklist.SV$CHR <- paste('chr', variants.blacklist.SV$CHR, sep = '')
variants.greenlist.SNV$CHR <- paste('chr', variants.greenlist.SNV$CHR, sep = '')
variants.greenlist.SV$CHR <- paste('chr', variants.greenlist.SV$CHR, sep = '')

# Subset the BED info
variants.blacklist.SNV <- variants.blacklist.SNV[, c('CHR', 'BP', 'BP', 'SNP')]
#variants.blacklist.SV <- variants.blacklist.SV[, c('CHR', 'BP', 'BP', 'SNP')]
variants.greenlist.SNV <- variants.greenlist.SNV[, c('CHR', 'BP', 'BP', 'SNP')]
variants.greenlist.SV <- variants.greenlist.SV[, c('CHR', 'BP', 'BP', 'SNP')]

    # Save the BED file (END POS IS THE SAME AS THE START POSITION)
    write.table(variants.blacklist.SNV, file = paste(dir.output, "Annotated_Significant_Variants_Blacklist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
    #write.table(variants.blacklist.SV, file = paste(dir.output, "Annotated_Significant_Variants_Blacklist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
    write.table(variants.greenlist.SNV, file = paste(dir.output, "Annotated_Significant_Variants_Greenlist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
    write.table(variants.greenlist.SV, file = paste(dir.output, "Annotated_Significant_Variants_Greenlist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         

# Make a second set of BED files with a window around each variant

    # Make variable to hold secondary positions
    secondary.black.SNV <- variants.blacklist.SNV
    secondary.green.SNV <- variants.greenlist.SNV
    secondary.green.SV <- variants.greenlist.SV
    
    # Update the positions (+/- 500 bp)
    secondary.black.SNV$BP <- secondary.black.SNV$BP - 500 
    secondary.green.SNV$BP <- secondary.green.SNV$BP - 500
    secondary.green.SV$BP <- secondary.green.SV$BP - 500
    
    secondary.black.SNV$BP.1 <- secondary.black.SNV$BP.1 + 500 
    secondary.green.SNV$BP.1 <- secondary.green.SNV$BP.1 + 500
    secondary.green.SV$BP.1 <- secondary.green.SV$BP.1 + 500
    
    # Save the BED file (POSITIONS HAVE +/- 500 BP WINDOW)
    write.table(secondary.black.SNV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_ALL_Variants_Blacklist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
    write.table(secondary.green.SNV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_ALL_Variants_Greenlist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
    write.table(secondary.green.SV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_ALL_Variants_Greenlist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)   
    
# Make a second set of BED files with 1) only index variants and 2) a window around each variant
    
    # Make variable to hold index variants
    sig.index.variants <- sig.variants
    
    # Subset index variants
    sig.index.variants <- sig.index.variants[which(sig.index.variants$Index_Variant == 'Yes'), ]
    
    # Get BED info for index variants
    variants.blacklist.SNV <- variants.blacklist.SNV[which(variants.blacklist.SNV$SNP %in% sig.index.variants$SNP), ]
    #variants.blacklist.SV <- variants.blacklist.SV[which(variants.blacklist.SV$SNP %in% sig.index.variants$SNP), ]
    variants.greenlist.SNV <- variants.greenlist.SNV[which(variants.greenlist.SNV$SNP %in% sig.index.variants$SNP), ]
    variants.greenlist.SV <- variants.greenlist.SV[which(variants.greenlist.SV$SNP %in% sig.index.variants$SNP), ]
    
        # Save the BED file (END POS IS THE SAME AS THE START POSITION)
        write.table(variants.blacklist.SNV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_INDEX_Variants_Blacklist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
        #write.table(variants.blacklist.SV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_INDEX_Variants_Blacklist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
        write.table(variants.greenlist.SNV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_INDEX_Variants_Greenlist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
        write.table(variants.greenlist.SV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_INDEX_Variants_Greenlist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         

    # Add window (+/- 500 bp)
    variants.blacklist.SNV$BP <- variants.blacklist.SNV$BP - 500
    #variants.blacklist.SV$BP <- variants.blacklist.SV$BP - 500
    variants.greenlist.SNV$BP <- variants.greenlist.SNV$BP - 500
    variants.greenlist.SV$BP <- variants.greenlist.SV$BP - 500
    
    variants.blacklist.SNV$BP.1 <- variants.blacklist.SNV$BP.1 + 500
    #variants.blacklist.SV$BP.1 <- variants.blacklist.SV$BP.1 + 500
    variants.greenlist.SNV$BP.1 <- variants.greenlist.SNV$BP.1 + 500
    variants.greenlist.SV$BP.1 <- variants.greenlist.SV$BP.1 + 500
    
        # Save the BED file (POSITIONS HAVE +/- 500 BP WINDOW)
        write.table(variants.blacklist.SNV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_INDEX_Variants_Blacklist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
        #write.table(variants.blacklist.SV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_INDEX_Variants_Blacklist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
        write.table(variants.greenlist.SNV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_INDEX_Variants_Greenlist_SNVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         
        write.table(variants.greenlist.SV, file = paste(dir.output.EXTRA.BED, "Annotated_Significant_WINDOWED_INDEX_Variants_Greenlist_SVs.bed", sep=""), sep = "\t" , row.names = F, col.names = F, quote = F)         

        
        
        
        
# END SCRIPT ------------------------------------------------------------------
        
        
# Clean the environment
rm(list=ls())       
        

