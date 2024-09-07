# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/5_ORA/Run_ORA_Functions.R")

# Load libraries
library(DOSE)
library(clusterProfiler) # to run ORA
library(ggplot2) # for split dotplots
library(scales) # for modifying the ggplot colorbar
library(enrichplot) # for gseaplot2

# Load gene sets 
load(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/5_ORA/Gene_Set_Collections_for_ORA_GSEA.R')



# ORA TETrancripts Consensus WGCNA ------------------------------------------------------------------------------------------------------------


# directory for output files
my.output.dir <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/5_ORA/Results_ORA_Consensus_GEUVADIS_WGCNA/'

# Load WGCNA results 
results.all <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Network_Construction_Consensus/TETranscripts_GEUVADIS_Consensus_Module_vs_Gene_Key.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = NULL)
rownames(results.all) <- results.all$gene

# define the universe genes
universe.entire <- as.character(rownames(results.all))

    # Keep only genes (i.e. remove TEs/Repeats)
    universe.entire <- universe.entire[grepl('ENSG', universe.entire)]

# Extract module labels
gene.modules <- unique(results.all$module)

# Loop over each cluster and run ORA
for (ith_cluster in gene.modules) {
  
    # Define the ith genelist (genes in the ith module)
    ith_gene_list <- rownames(results.all[which(results.all$module == ith_cluster), ])
    
    # using GO BP Collection
    ORA.GO <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'GO/', universe.list = universe.entire, my.genelist = ith_gene_list, my.gs.collection = GO_BP_Geneset, gs.label = 'GO_BP', condition_label = paste('Module_', ith_cluster, sep = ''))

    # using Hallmark
    ORA.Hallmark <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'Hallmark/', universe.list = universe.entire, my.genelist = ith_gene_list, my.gs.collection = Hallmark_Geneset, gs.label = 'Hallmark', condition_label = paste('Module_', ith_cluster, sep = ''))

     # using Reactome
    ORA.Reactome <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'Reactome/', universe.list = universe.entire, my.genelist = ith_gene_list, my.gs.collection = Reactome_Geneset, gs.label = 'Reactome', condition_label = paste('Module_', ith_cluster, sep = ''))

    # using miRDB
    ORA.miRDB <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'miRDB/', universe.list = universe.entire, my.genelist = ith_gene_list, my.gs.collection = miRDB_Geneset, gs.label = 'miRDB', condition_label = paste('Module_', ith_cluster, sep = ''))

    # using TFT GTRD
    ORA.TFT <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'TFT_GTRD/', universe.list = universe.entire, my.genelist = ith_gene_list, my.gs.collection = GTRD_Geneset, gs.label = 'TFT_GTRD', condition_label = paste('Module_', ith_cluster, sep = ''))

    # using CRISPR screen regulators
    ORA.regulators <- Run_ORA(output.dir = my.output.dir, output_subfolder = 'L1_CRISPR_Screen_Regulators/', universe.list = universe.entire, my.genelist = ith_gene_list, my.gs.collection = L1_regulator.genesets, gs.label = 'L1_Regulators', condition_label = paste('Module_', ith_cluster, sep = ''))

}



# Save Session Info  ------------------------------------------------------------------------------------------------------------


# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/5_ORA/Session_Info/'
    
    
sink(file = paste(dir.session_info, "Session_Info_ORA_Consensus_WGCNA.txt", sep =""))
sessionInfo()
sink()      
    

    
# Clean the environment
rm(list=ls())
    
    
