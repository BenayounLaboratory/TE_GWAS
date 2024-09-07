# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Run_ORA_Functions.R")

# Load libraries
library(msigdbr)
library(DOSE)
library(clusterProfiler) # to run ORA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2

# RETRIEVE CURATED GENE SETS  

    # Retrieve MSigDBR GO BP gene sets
    GO_BP_Geneset <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = 'GO:BP') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR GO CC gene sets
    GO_CC_Geneset <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = 'GO:CC') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR GO MF gene sets
    GO_MF_Geneset <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = 'GO:MF') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR HALLMARK gene sets
    Hallmark_Geneset <- msigdbr(species = "Homo sapiens", category = 'H') %>% 
    dplyr::select(gs_name, gene_symbol)
      
    # Retrieve MSigDBR REACTOME gene sets
    Reactome_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:REACTOME') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR Biocarta gene sets
    Biocarta_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:BIOCARTA') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR PID gene sets
    PID_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:PID') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR Wikipathways gene sets
    Wikipathways_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:WIKIPATHWAYS') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR KEGG gene sets
    KEGG_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CP:KEGG') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR miRDB gene set
    miRDB_Geneset <- msigdbr(species = "Homo sapiens", category = 'C3', subcategory = 'MIR:MIRDB') %>% 
    dplyr::select(gs_name, gene_symbol)
          
    # Retrieve MSigDBR TFT GTRD gene set
    GTRD_Geneset <- msigdbr(species = "Homo sapiens", category = 'C3', subcategory = 'TFT:GTRD') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve MSigDBR CGP (chemical and genetic perturbations) gene set
    CGP_Geneset <- msigdbr(species = "Homo sapiens", category = 'C2', subcategory = 'CGP') %>% 
    dplyr::select(gs_name, gene_symbol)
    
    # Retrieve Interpro Domain gene set
    Interpro_Geneset <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/External_Genesets/Cleaned_Up_Interpro_Domain_Geneset.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
    Interpro_Geneset <- Interpro_Geneset[, c('V4', 'V5')]
    colnames(Interpro_Geneset) <- c('gs_name', 'gene_symbol')
    
    

# SNVs Greenlist + Blacklist ------------------------------------------------------------------------------------------------------------


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Overrepresentation_Analysis/SNVs_Greenlist_and_Blacklist/'

# Load genelist for SNVs
genes.greenlist.SNVs <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_Significant_Variants_Greenlist_SNVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
genes.greenlist.SNVs <- genes.greenlist.SNVs$V1

genes.blacklist.SNVs <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_Significant_Variants_Blacklist_SNVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
genes.blacklist.SNVs <- genes.blacklist.SNVs$V1

# Combine black and green list genes, keeping only unique genes
my.sig.genes <- c(genes.greenlist.SNVs, genes.blacklist.SNVs)
my.sig.genes <- unique(my.sig.genes)

# Load genelist for the universe
my.universe <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_BACKGROUND_SNVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
my.universe <- my.universe$V1
         
# Run ORA
    
    # using GO BP Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using GO CC Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = GO_CC_Geneset, gs.label = 'GO_CC', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using GO MF Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = GO_MF_Geneset, gs.label = 'GO_MF', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using Hallmark
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using Reactome
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Reactome/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = Reactome_Geneset, gs.label = 'Reactome', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using Biocarta
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Biocarta/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = Biocarta_Geneset, gs.label = 'Biocarta', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using PID
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'PID/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = PID_Geneset, gs.label = 'PID', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using Wikipathways
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Wikipathways/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = Wikipathways_Geneset, gs.label = 'Wikipathways', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))
    
    # using KEGG
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'KEGG/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = KEGG_Geneset, gs.label = 'KEGG', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))
    
    # using miRDB
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'miRDB/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using TFT GTRD
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = GTRD_Geneset, gs.label = 'TFT_GTRD', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using CGP
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'CGP/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = CGP_Geneset, gs.label = 'CGP', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))

    # using Interpro Domains
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Interpro_Domains/', universe.list = my.universe, my.genelist = my.sig.genes, my.gs.collection = Interpro_Geneset, gs.label = 'Interpro_Domains', condition_label = paste('Black_and_Greenlist_SNVs', sep = ''))


# SNVs Greenlist ------------------------------------------------------------------------------------------------------------


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Overrepresentation_Analysis/SNVs_Greenlist/'

# Load genelist for greenlist SNVs
genes.greenlist.SNVs <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_Significant_Variants_Greenlist_SNVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
genes.greenlist.SNVs <- genes.greenlist.SNVs$V1

# Load genelist for the universe
my.universe <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_BACKGROUND_SNVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
my.universe <- my.universe$V1
         
# Run ORA
    
    # using GO BP Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using GO CC Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = GO_CC_Geneset, gs.label = 'GO_CC', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using GO MF Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = GO_MF_Geneset, gs.label = 'GO_MF', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using Hallmark
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using Reactome
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Reactome/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = Reactome_Geneset, gs.label = 'Reactome', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using Biocarta
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Biocarta/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = Biocarta_Geneset, gs.label = 'Biocarta', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using PID
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'PID/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = PID_Geneset, gs.label = 'PID', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using Wikipathways
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Wikipathways/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = Wikipathways_Geneset, gs.label = 'Wikipathways', condition_label = paste('Greenlist_SNVs', sep = ''))
    
    # using KEGG
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'KEGG/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = KEGG_Geneset, gs.label = 'KEGG', condition_label = paste('Greenlist_SNVs', sep = ''))
    
    # using miRDB
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'miRDB/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using TFT GTRD
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = GTRD_Geneset, gs.label = 'TFT_GTRD', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using CGP
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'CGP/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = CGP_Geneset, gs.label = 'CGP', condition_label = paste('Greenlist_SNVs', sep = ''))

    # using Interpro Domains
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Interpro_Domains/', universe.list = my.universe, my.genelist = genes.greenlist.SNVs, my.gs.collection = Interpro_Geneset, gs.label = 'Interpro_Domains', condition_label = paste('Greenlist_SNVs', sep = ''))


# SVs Greenlist ------------------------------------------------------------------------------------------------------------


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Overrepresentation_Analysis/SVs_Greenlist/'

# Load genelist for greenlist SNVs
genes.greenlist.SVs <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_Significant_Variants_Greenlist_SVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
genes.greenlist.SVs <- genes.greenlist.SVs$V1

# Load genelist for the universe
my.universe <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Annotated_Significant_Variants/NEARBY_GENE_LIST_BACKGROUND_SVs.txt", header = F, stringsAsFactors = F, sep = '\t', row.names = NULL)
my.universe <- my.universe$V1
         
# Run ORA
    
    # using GO BP Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = paste('Greenlist_SVs', sep = ''))

    # using GO CC Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = GO_CC_Geneset, gs.label = 'GO_CC', condition_label = paste('Greenlist_SVs', sep = ''))

    # using GO MF Collection
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = GO_MF_Geneset, gs.label = 'GO_MF', condition_label = paste('Greenlist_SVs', sep = ''))

    # using Hallmark
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = paste('Greenlist_SVs', sep = ''))

    # using Reactome
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Reactome/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = Reactome_Geneset, gs.label = 'Reactome', condition_label = paste('Greenlist_SVs', sep = ''))

    # using Biocarta
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Biocarta/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = Biocarta_Geneset, gs.label = 'Biocarta', condition_label = paste('Greenlist_SVs', sep = ''))

    # using PID
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'PID/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = PID_Geneset, gs.label = 'PID', condition_label = paste('Greenlist_SVs', sep = ''))

    # using Wikipathways
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Wikipathways/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = Wikipathways_Geneset, gs.label = 'Wikipathways', condition_label = paste('Greenlist_SVs', sep = ''))
    
    # using KEGG
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'KEGG/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = KEGG_Geneset, gs.label = 'KEGG', condition_label = paste('Greenlist_SVs', sep = ''))
    
    # using miRDB
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'miRDB/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = paste('Greenlist_SVs', sep = ''))

    # using TFT GTRD
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = GTRD_Geneset, gs.label = 'TFT_GTRD', condition_label = paste('Greenlist_SVs', sep = ''))

    # using CGP
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'CGP/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = CGP_Geneset, gs.label = 'CGP', condition_label = paste('Greenlist_SVs', sep = ''))

    # using Interpro Domains
    Run_ORA(output.dir = my.output.dir, output_subfolder = 'Interpro_Domains/', universe.list = my.universe, my.genelist = genes.greenlist.SVs, my.gs.collection = Interpro_Geneset, gs.label = 'Interpro_Domains', condition_label = paste('Greenlist_SVs', sep = ''))


# Save Session Info  ------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/5_Annotation_Enrichments/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_ORA_Significant_Variants.txt", sep =""))
sessionInfo()
sink()      
    

    
# Clean the environment
rm(list=ls())
    
    
