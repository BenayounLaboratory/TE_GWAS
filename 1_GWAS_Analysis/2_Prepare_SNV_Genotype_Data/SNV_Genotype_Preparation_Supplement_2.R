# Set strings as factors
options(stringsAsFactors = F)

# Load libraries

# Specify output dir
dir.genotypes.AFR <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/'



# MAKE VARIANT MAPS AND TXT/BED FILES FOR SV DATA -------------------------------------------------------------------------------------------------


# GENERATE POSITION/MAPPING FILES

# Load initial data from plink BIM file
snp_info.AFR <- read.csv('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Structural_Variant_Genotypes/Filter_2/Plink_AFR_Genotypes/Final_SVs_AFR.bim', header = F, stringsAsFactors = F, sep = '\t')

# Rename columns. !!!!!!!!!!!!!!!!!!! NOTE: If -keep-allele-order is used in plink, REF allele is automatically assigned to A2 (last column in BIM file). NOTE: The REF allele may or may not be the Major allele.
colnames(snp_info.AFR) <- c('chr', 'ID', 'Morgan_Pos', 'pos', 'ALT', 'REF')

# Check if 'ID' labels are unique
length(unique(snp_info.AFR$ID)) == nrow(snp_info.AFR) # They are all unique. The code below may need to be modified if the pipeline is changed and IDs end up not being unique.

# Add snpid column. Format: row#_chr_pos. Important since some snps can fall in the same position (though they don't in this analysis)
snp_info.AFR$snpid <- paste(rownames(snp_info.AFR), snp_info.AFR$chr, snp_info.AFR$pos, sep = '_')

# Rearrange columns, removing "Morgan position"
snp_info.AFR <- snp_info.AFR[, c('snpid', 'chr', 'pos', 'ID', 'REF', 'ALT')]

    # Save table
    write.table(snp_info.AFR, file = paste(dir.genotypes.AFR, "RESOURCE_Variant_Map_AFR_SVs_only", ".txt", sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    # Save a list of SV IDs
    write.table(snp_info.AFR$ID, file = paste(dir.genotypes.AFR, "RESOURCE_IDs_SVs_only", ".txt", sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)
    
    
    
# GENERATE BED FILE
    
# Make a variable to hold the BED info    
my.snp.bed <- snp_info.AFR[, c('chr', 'pos', 'pos', 'ID')]

# Update the chromosome labels to include 'chr' (needed by pipelines like GREAT)
my.snp.bed$chr <- paste('chr', my.snp.bed$chr, sep = '')
      
# Save the BED file (END POS IS THE SAME AS THE START POSITION)
write.table(my.snp.bed, file = paste(dir.genotypes.AFR, "RESOURCE_All_Background_Positions_SVs_only", ".bed", sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)





# MAKE VARIANT MAPS AND TXT/BED FILES FOR SNV DATA -------------------------------------------------------------------------------------------------


# GENERATE POSITION/MAPPING FILES

# Load initial data from plink BIM file
snp_info.AFR <- read.csv('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/SNV_Genotypes/Filtered_2/Plink_AFR_Genotypes/Final_SNVs_AFR.bim', header = F, stringsAsFactors = F, sep = '\t')

# Rename columns. !!!!!!!!!!!!!!!!!!! NOTE: If -keep-allele-order is used in plink, REF allele is automatically assigned to A2 (last column in BIM file). NOTE: The REF allele may or may not be the Major allele.
colnames(snp_info.AFR) <- c('chr', 'RsID', 'Morgan_Pos', 'pos', 'ALT', 'REF')

# Check if 'RsID' labels are unique
length(unique(snp_info.AFR$RsID)) == nrow(snp_info.AFR) # They are all unique. The code below may need to be modified if the pipeline is changed and IDs end up not being unique.

# Add snpid column. Format: row#_chr_pos. Important since some snps can fall in the same position (though they don't in this analysis)
snp_info.AFR$snpid <- paste(rownames(snp_info.AFR), snp_info.AFR$chr, snp_info.AFR$pos, sep = '_')

# Rearrange columns, removing "Morgan position"
snp_info.AFR <- snp_info.AFR[, c('snpid', 'chr', 'pos', 'RsID', 'REF', 'ALT')]

    # Save table
    write.table(snp_info.AFR, file = paste(dir.genotypes.AFR, "RESOURCE_Variant_Map_AFR_SNVs_only", ".txt", sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    # Save a list of SNV RsIDs
    write.table(snp_info.AFR$RsID, file = paste(dir.genotypes.AFR, "RESOURCE_IDs_SNVs_only", ".txt", sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)


    
    
    
# GENERATE BED FILE
    
# Make a variable to hold the BED info    
my.snp.bed <- snp_info.AFR[, c('chr', 'pos', 'pos', 'RsID')]

# Update the chromosome labels to include 'chr' (needed by pipelines like GREAT)
my.snp.bed$chr <- paste('chr', my.snp.bed$chr, sep = '')
      
# Save the BED file (END POS IS THE SAME AS THE START POSITION)
write.table(my.snp.bed, file = paste(dir.genotypes.AFR, "RESOURCE_All_Background_Positions_SNVs_only", ".bed", sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)





# MAKE VARIANT MAPS AND TXT/BED FILES FOR SNV+SV DATA -------------------------------------------------------------------------------------------------


# GENERATE POSITION/MAPPING FILES

# Load initial data from plink BIM file
snp_info.AFR <- read.csv('/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/Final_SNVs_SVs_AFR.bim', header = F, stringsAsFactors = F, sep = '\t')

# Rename columns. !!!!!!!!!!!!!!!!!!! NOTE: If -keep-allele-order is used in plink, REF allele is automatically assigned to A2 (last column in BIM file). NOTE: The REF allele may or may not be the Major allele.
colnames(snp_info.AFR) <- c('chr', 'RsID', 'Morgan_Pos', 'pos', 'ALT', 'REF')

# Check if 'RsID' labels are unique
length(unique(snp_info.AFR$RsID)) == nrow(snp_info.AFR) # They are all unique. The code below may need to be modified if the pipeline is changed and IDs end up not being unique.

# Add snpid column. Format: row#_chr_pos. Important since some snps can fall in the same position (though they don't in this analysis)
snp_info.AFR$snpid <- paste(rownames(snp_info.AFR), snp_info.AFR$chr, snp_info.AFR$pos, sep = '_')

# Rearrange columns, removing "Morgan position"
snp_info.AFR <- snp_info.AFR[, c('snpid', 'chr', 'pos', 'RsID', 'REF', 'ALT')]

    # Save table
    write.table(snp_info.AFR, file = paste(dir.genotypes.AFR, "RESOURCE_Variant_Map_AFR_SNVs_and_SVs", ".txt", sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    # Save a list of SNV/SV IDs
    write.table(snp_info.AFR$RsID, file = paste(dir.genotypes.AFR, "RESOURCE_IDs_SNVs_and_SVs", ".txt", sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)
    

    
    
        
# GENERATE SNV+SV BED FILE
    
# Make a variable to hold the BED info    
my.snp.bed <- snp_info.AFR[, c('chr', 'pos', 'pos', 'RsID')]

# Update the chromosome labels to include 'chr' (needed by pipelines like GREAT)
my.snp.bed$chr <- paste('chr', my.snp.bed$chr, sep = '')
      
# Save the BED file (END POS IS THE SAME AS THE START POSITION)
write.table(my.snp.bed, file = paste(dir.genotypes.AFR, "RESOURCE_All_Background_Positions_SNVs_and_SVs", ".bed", sep =""), sep = "\t" , row.names = F, col.names = F, quote = F)





# END -------------------------------------------------------------------------------------------------


# Remove unneeded variables
rm(list=ls()) 
      
