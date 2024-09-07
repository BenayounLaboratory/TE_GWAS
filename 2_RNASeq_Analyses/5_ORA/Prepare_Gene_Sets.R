# Set strings as factors
options(stringsAsFactors = F)

# Load libraries
library(msigdbr) # Make sure this is version 7.5.1; this contains MSigDB genesets with genes in ENSEMBL format
library(stringr) # for capitalization changes
library(splitstackshape)

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/5_ORA/'             
            
            

# GENERATE ENSEMBL GENE SETS ----------------------------------------------------------------------------------------------------



# RETRIEVE CURATED GENE SETS  

# Retrieve MSigDBR HALLMARK gene sets
Hallmark_Geneset <- msigdbr(species = "Homo sapiens", category = 'H') %>% 
dplyr::select(gs_name, human_ensembl_gene)
  
# Retrieve MSigDBR REACTOME gene sets
Reactome_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:REACTOME') %>% 
dplyr::select(gs_name, human_ensembl_gene)

# Retrieve MSigDBR GO BP gene sets
GO_BP_Geneset <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = 'GO:BP') %>% 
dplyr::select(gs_name, human_ensembl_gene)

# Retrieve MSigDBR miRDB gene set
miRDB_Geneset <- msigdbr(species = "Homo sapiens", category = 'C3', subcategory = 'MIR:MIRDB') %>% 
dplyr::select(gs_name, human_ensembl_gene)
      
# Retrieve MSigDBR TFT GTRD gene set
GTRD_Geneset <- msigdbr(species = "Homo sapiens", category = 'C3', subcategory = 'TFT:GTRD') %>% 
dplyr::select(gs_name, human_ensembl_gene)
         




# CLEANUP GENE SET NAMES 

# remove redundant parts of labels
Hallmark_Geneset$gs_name <- gsub("HALLMARK_", "", Hallmark_Geneset$gs_name)
Reactome_Geneset$gs_name <- gsub("REACTOME_", "", Reactome_Geneset$gs_name)
GO_BP_Geneset$gs_name <- gsub("GOBP_", "", GO_BP_Geneset$gs_name)

# write function to lowercase all letters, then uppercase the first
sentence_casing <- function(input_names){
  
  # lowercase all letters
  input_lowercased <- tolower(input_names)
  
  # uppercase the first letter
  input_sentence_case <- str_to_sentence(input_lowercased)
  
  # output the final string
  return(input_sentence_case)
  
}

# change cases of gene set names, so only the first letter is capitalized
Hallmark_Geneset$gs_name <- sentence_casing(input_names = Hallmark_Geneset$gs_name)
Reactome_Geneset$gs_name <- sentence_casing(input_names = Reactome_Geneset$gs_name)
GO_BP_Geneset$gs_name <- sentence_casing(input_names = GO_BP_Geneset$gs_name)
GTRD_Geneset$gs_name <- sentence_casing(input_names = GTRD_Geneset$gs_name)

# Update underscores to either blank spaces or dashes
Hallmark_Geneset$gs_name <- gsub("_", " ", Hallmark_Geneset$gs_name)
Reactome_Geneset$gs_name <- gsub("_", " ", Reactome_Geneset$gs_name)
GO_BP_Geneset$gs_name <- gsub("_", " ", GO_BP_Geneset$gs_name)
miRDB_Geneset$gs_name <- gsub("_", "-", miRDB_Geneset$gs_name)
GTRD_Geneset$gs_name <- gsub("_", " ", GTRD_Geneset$gs_name)


  


# MAKE TRANSPOSON GENE SETS (BY CLASS AND FAMILY, USING TETRANSCRIPTS DESIGNATIONS)

# Load TETranscripts raw counts file (which has all the annotated transposons)  
my.counts <- read.csv(file = "/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/Quantifications_from_Bravo_et_al_2024/TETranscripts_Counts_EUR/fastp_ERR188021Aligned.sortedByCoord.out.bam.cntTable", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
  
# Define rows that are TEs
TE.rows <- which(!grepl('ENSG', rownames(my.counts)))

# Make dataframe to hold TE classifications
my.TE.classifications <- data.frame(row.names = rownames(my.counts[TE.rows, , drop = FALSE]),
                                    Full_label = rownames(my.counts[TE.rows, , drop = FALSE]),
                                    To_split = rownames(my.counts[TE.rows, , drop = FALSE]))
  
# Split TE labels into subfamily, family, and class. 
my.TE.classifications <- as.data.frame(splitstackshape::cSplit(my.TE.classifications, splitCols = "To_split", sep = ":", direction = "wide", fixed = TRUE, drop = TRUE, stripWhite = FALSE, makeEqual = TRUE))

# Update the column names
names(my.TE.classifications)[2:4] <- c('Subfamily', 'Family', 'Class')   

# Append "subfamilies" to the Family/Class names, since these names will be used as the gene set names
my.TE.classifications$Class <- paste(my.TE.classifications$Class, ' subfamilies', sep = '')
my.TE.classifications$Family <- paste(my.TE.classifications$Family, ' subfamilies', sep = '')
  
# Make gene set (ready for GSEA)
Repeat_class.gs <- data.frame(gs_name = my.TE.classifications$Class, gene = my.TE.classifications$Full_label)
Repeat_family.gs <- data.frame(gs_name = my.TE.classifications$Family, gene = my.TE.classifications$Full_label)
           
            
        


# MAKE Alu GENE SETS (BY Alu AGE, USING TETRANSCRIPTS DESIGNATIONS)
  
# Load TETranscripts raw counts file (which has all the annotated transposons)  
my.counts <- read.csv(file = "/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/0_All_Sample_Metadata_and_External_Resources/Quantifications_from_Bravo_et_al_2024/TETranscripts_Counts_EUR/fastp_ERR188021Aligned.sortedByCoord.out.bam.cntTable", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
  
# Define rows that are Alu
L1.rows <- which(grepl(':Alu:', rownames(my.counts)))

# Make dataframe to hold TE classifications
my.Alu.classifications <- data.frame(row.names = rownames(my.counts[L1.rows, , drop = FALSE]),
                                    Full_label = rownames(my.counts[L1.rows, , drop = FALSE]),
                                    To_split = rownames(my.counts[L1.rows, , drop = FALSE]))


# Segregate by age
AluJ_subfamilies <- my.Alu.classifications[which(grepl('AluJ', rownames(my.Alu.classifications))), ]
AluS_subfamilies <- my.Alu.classifications[which(grepl('AluS', rownames(my.Alu.classifications))), ]
AluY_subfamilies <- my.Alu.classifications[which(grepl('AluY', rownames(my.Alu.classifications))), ]

# Add age columns
AluJ_subfamilies$Age <- 'AluJ subfamilies'
AluS_subfamilies$Age <- 'AluS subfamilies'
AluY_subfamilies$Age <- 'AluY subfamilies'

# Rejoin Alu subfamilies
my.Alu.classifications <- rbind(AluJ_subfamilies, AluS_subfamilies, AluY_subfamilies)
  
# Make gene set (ready for GSEA)
Alu_by_age.gs <- data.frame(gs_name = my.Alu.classifications$Age, gene = my.Alu.classifications$Full_label)





# MAKE L1 CRISPR SCREEN REGULATOR GENESETS

# Load lists
regulator.expr <- read.csv(file = "/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Li_2024_L1_Expression_Regulators/Li_2024_L1_Expression_Regulators_Both_Activators_Suppressors.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
regulator.transpose.primary <- read.csv(file = "/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Liu_2018_Transposition_Regulators/L1_regulators_Primary_Screen.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')
regulator.transpose.optimized <- read.csv(file = "/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Liu_2018_Transposition_Regulators/L1_regulators_L1opt_Screen.txt", header = TRUE, row.names = 1, stringsAsFactors = FALSE, sep = '\t')

    # Remove genes with NA
    regulator.expr <- na.omit(regulator.expr)
    regulator.transpose.primary <- na.omit(regulator.transpose.primary)
    regulator.transpose.optimized <- na.omit(regulator.transpose.optimized)
    
    # Add Ensembl column
    regulator.transpose.primary$Ensembl <- rownames(regulator.transpose.primary)
    regulator.transpose.optimized$Ensembl <- rownames(regulator.transpose.optimized)

# Keep only the columns necessary for the gene set
regulator.expr <- regulator.expr[, c('Effect', 'gProfiler_Ensembl')]
regulator.transpose.primary <- regulator.transpose.primary[, c('Cell_Effect', 'Ensembl')]
regulator.transpose.optimized <- regulator.transpose.optimized[, c('Sig', 'Ensembl')]

# Use common column names
colnames(regulator.expr) <- c('gs_name', 'gene')
colnames(regulator.transpose.primary) <- c('gs_name', 'gene')
colnames(regulator.transpose.optimized) <- c('gs_name', 'gene')

# Update geneset names
regulator.expr[regulator.expr == 'Activator'] <- 'L1 Expression Activators'
regulator.expr[regulator.expr == 'Suppressor'] <- 'L1 Expression Suppresors'

regulator.transpose.primary[regulator.transpose.primary == 'Both_Act'] <- 'L1 Transposition Activators'
regulator.transpose.primary[regulator.transpose.primary == 'Both_Supp'] <- 'L1 Transposition Suppressors'
regulator.transpose.primary[regulator.transpose.primary == 'HeLa_Act'] <- 'L1 Transposition Regulators (varying effects)'
regulator.transpose.primary[regulator.transpose.primary == 'HeLa_Supp'] <- 'L1 Transposition Regulators (varying effects)'
regulator.transpose.primary[regulator.transpose.primary == 'K562_Act'] <- 'L1 Transposition Regulators (varying effects)'
regulator.transpose.primary[regulator.transpose.primary == 'K562_Supp'] <- 'L1 Transposition Regulators (varying effects)'
regulator.transpose.primary[regulator.transpose.primary == 'KA_HS'] <- 'L1 Transposition Regulators (varying effects)'
regulator.transpose.primary[regulator.transpose.primary == 'KS_HA'] <- 'L1 Transposition Regulators (varying effects)'

regulator.transpose.optimized[regulator.transpose.optimized == 'opt-L1_Act'] <- 'Codon-Optimized-L1 Transposition Activators'
regulator.transpose.optimized[regulator.transpose.optimized == 'opt-L1_Supp'] <- 'Codon-Optimized-L1 Transposition Suppressors'
 
# Combine genesets
L1_regulator.genesets <- rbind(regulator.expr, regulator.transpose.primary, regulator.transpose.optimized)

    # Update rownames
    rownames(L1_regulator.genesets) <- 1:nrow(L1_regulator.genesets)
     
    
    
    
    
# Save gene sets
save(Hallmark_Geneset,
     Reactome_Geneset,
     GO_BP_Geneset,
     miRDB_Geneset,
     GTRD_Geneset,
     Repeat_class.gs,
     Repeat_family.gs,
     Alu_by_age.gs,
     L1_regulator.genesets,
     file = paste(dir.output, "Gene_Set_Collections_for_ORA_GSEA.R", sep = ''))
          




# SESSION INFO ----------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/5_ORA/Session_Info/'
    
sink(file = paste(dir.session_info, "Session_Info_Prepare_Gene_Sets.txt", sep =""))
sessionInfo()
sink()   

      
# Clean the environment
rm(list=ls())      