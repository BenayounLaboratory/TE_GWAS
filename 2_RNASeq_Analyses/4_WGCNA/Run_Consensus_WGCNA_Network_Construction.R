# Set strings as factors
options(stringsAsFactors = F)

# Load functions associated with this script.
source("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/WGCNA_Functions.R") 

# Load libraries
library(WGCNA) # for constructing the gene network
library(ggplot2) # for making plots
library(cowplot) # for ggplot image grids?
library(gridExtra) # for ggplot image grids?
library(pheatmap) # for heatmaps
library(biomaRt) # for mapping Ensembl to other IDs

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Network_Construction_Consensus/'

# RESOURCES:
# Tutorial from https://alexslemonade.github.io/refinebio-examples/04-advanced-topics/network-analysis_rnaseq_01_wgcna.html
# Tutorial from WGCNA website tutorials
# Theory and consideration: https://edu.isb-sib.ch/pluginfile.php/158/course/section/65/_01_SIB2016_wgcna.pdf



# RUN WGCNA (TETranscripts data) --------------------------------------------


# FIND SOFT THRESHOLD PARAMETER
  
# Define output directory
setwd(dir.output)

# Load expression data.
VST.expression.EUR <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/All_counts_EUR_358_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_EBV.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)
VST.expression.YRI <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/2_Prepare_RNASeq/Processed_counts/All_counts_YRI_86_filtered_VST_BatchesRemoved_lab_ancestry_PCs_sex_EBV.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Transpose the data to match WGCNA input requirements
    VST.expression.EUR <- t(VST.expression.EUR)
    VST.expression.YRI <- t(VST.expression.YRI)
    
# Put the expression data into a multiset for consensus network generation. 
my.multiExpr = vector(mode = "list", length = 2)

    
    my.multiExpr[[1]] = list(data = as.data.frame(VST.expression.EUR))
    names(my.multiExpr[[1]]$data) = colnames(VST.expression.EUR)
    rownames(my.multiExpr[[1]]$data) = rownames(VST.expression.EUR)
    
    my.multiExpr[[2]] = list(data = as.data.frame(VST.expression.YRI))
    names(my.multiExpr[[2]]$data) = colnames(VST.expression.YRI)
    rownames(my.multiExpr[[2]]$data) = rownames(VST.expression.YRI)
    
    # Check that the data has the correct format for many functions operating on multiple sets:
    checkSets(my.multiExpr)
    
# Adjacency thresholding parameter
my.sft.EUR <- pickSoftThreshold(VST.expression.EUR,
                                 dataIsExpr = TRUE,
                                 corFnc = bicor,
                                 corOptions = list(use = 'p', maxPOutliers = 0.05),
                                 networkType = "signed",
                                 powerVector = seq(1, 20, by = 1),
                                 verbose = 5,
                                 blockSize = 20000
                                )

my.sft.YRI <- pickSoftThreshold(VST.expression.YRI,
                                 dataIsExpr = TRUE,
                                 corFnc = bicor,
                                 corOptions = list(use = 'p', maxPOutliers = 0.05),
                                 networkType = "signed",
                                 powerVector = seq(1, 20, by = 1),
                                 verbose = 5,
                                 blockSize = 20000
                                )


# PLOT SCALE INDEPENDENCE (EUR)

    # Extract data from sft object in df format
    my.sft_df <- data.frame(my.sft.EUR$fitIndices)
    
    # calculate model fit
    my.sft_df$model_fit <- -sign(my.sft_df$slope) * my.sft_df$SFT.R.sq

    # Start plot
    p1 <- ggplot(my.sft_df, aes(x = Power, y = model_fit, label = Power)) +

    # Plot the points
    geom_point() +
      
    # We'll put the Power labels slightly above the data points
    geom_text(nudge_y = 0.1) +
      
    # WGCNA authors recommend an R^2 greater than 0.80. You want a value that is sufficiently high, but not too high
    #geom_hline(yintercept = 0.80, col = "red") +
    geom_hline(yintercept = 0.90, col = "red") +
      
    # Set the plot xlims
    xlim(0, 20) +

    # Just in case our values are low, we want to make sure we can still see the 0.80 level
    ylim(c(min(my.sft_df$model_fit), 1.10)) +
      
    # We can add more sensible labels for our axis
    xlab("Soft Threshold (power)") +
    ylab("Scale Free Topology Model Fit, signed R^2") +
    ggtitle("Scale independence") +
      
    # This adds some nicer aesthetics to our plot
    theme_classic()
    
# PLOT MEAN CONNECTIVITY
    
    # Start Plot
    p2 <- ggplot(my.sft_df, aes(x = Power, y = mean.k., label = Power)) +
      
    # Plot the points
    geom_point() +
      
    # We'll put the Power labels slightly above the data points
    geom_text(nudge_y = 100) +
      
    # Set the plot xlims
    xlim(0, 20) +
      
    # Set the plot ylims
    ylim(0, 1000) +
      
    # arbitrary mean connectivity threshold
    geom_hline(yintercept = 100, col = "red") +
      
     # We can add more sensible labels for our axis
    xlab("Soft Threshold (power)") +
    ylab("Mean Connectivity") +
    ggtitle("Mean Connectivity") +
      
    # This adds some nicer aesthetics to our plot
    theme_classic()

    # Save file
    pdf(paste(dir.output, "TETranscripts_GEUVADIS_EUR_Power_vs_Scale_Independence_&_Mean_Connectivity", ".pdf", sep=""), width = 15, height = 15)
    plot_grid(p1, p2, nrow = 2, ncol = 2)
    dev.off()
    
    
# PLOT SCALE INDEPENDENCE (YRI)

    # Extract data from sft object in df format
    my.sft_df <- data.frame(my.sft.YRI$fitIndices)
    
    # calculate model fit
    my.sft_df$model_fit <- -sign(my.sft_df$slope) * my.sft_df$SFT.R.sq

    # Start plot
    p1 <- ggplot(my.sft_df, aes(x = Power, y = model_fit, label = Power)) +

    # Plot the points
    geom_point() +
      
    # We'll put the Power labels slightly above the data points
    geom_text(nudge_y = 0.1) +
      
    # WGCNA authors recommend an R^2 greater than 0.80. You want a value that is sufficiently high, but not too high
    #geom_hline(yintercept = 0.80, col = "red") +
    geom_hline(yintercept = 0.90, col = "red") +
      
    # Set the plot xlims
    xlim(0, 20) +

    # Just in case our values are low, we want to make sure we can still see the 0.80 level
    ylim(c(min(my.sft_df$model_fit), 1.10)) +
      
    # We can add more sensible labels for our axis
    xlab("Soft Threshold (power)") +
    ylab("Scale Free Topology Model Fit, signed R^2") +
    ggtitle("Scale independence") +
      
    # This adds some nicer aesthetics to our plot
    theme_classic()
    
# PLOT MEAN CONNECTIVITY
    
    # Start Plot
    p2 <- ggplot(my.sft_df, aes(x = Power, y = mean.k., label = Power)) +
      
    # Plot the points
    geom_point() +
      
    # We'll put the Power labels slightly above the data points
    geom_text(nudge_y = 100) +
      
    # Set the plot xlims
    xlim(0, 20) +
      
    # Set the plot ylims
    ylim(0, 1000) +
      
    # arbitrary mean connectivity threshold
    geom_hline(yintercept = 100, col = "red") +
      
     # We can add more sensible labels for our axis
    xlab("Soft Threshold (power)") +
    ylab("Mean Connectivity") +
    ggtitle("Mean Connectivity") +
      
    # This adds some nicer aesthetics to our plot
    theme_classic()

    # Save file
    pdf(paste(dir.output, "TETranscripts_GEUVADIS_YRI_Power_vs_Scale_Independence_&_Mean_Connectivity", ".pdf", sep=""), width = 15, height = 15)
    plot_grid(p1, p2, nrow = 2, ncol = 2)
    dev.off()

# Estimate Power (if you're not sure from the plots.)
power_estimate.EUR = my.sft.EUR$powerEstimate 
power_estimate.YRI = my.sft.YRI$powerEstimate 

# Set Power
softPower.EUR = 12
softPower.YRI = 12
softPower.consensus = 12

        
        
        
        
            
# AUTOMATIC NETWORK CONSTRUCTION         

# Find coexpression modules
my.net.consensus <- blockwiseConsensusModules(my.multiExpr, 
                                              corType = 'bicor',
                                              power = softPower.consensus, # soft threshold for network construction
                                              networkType = "signed",
                                              maxPOutliers = 0.05,
                                              maxBlockSize = 50000, # What size chunks (how many genes) the calculations should be run in. 30000 for 32GB computer
                                              TOMType = "signed", # topological overlap matrix
                                              saveTOMs = TRUE, 
                                              saveTOMFileBase = "TETranscripts_GEUVADIS_CONSENSUS_TOM",
                                              mergeCutHeight = 0.25, # for joining nearby modules. Modules with correlation coefficient 1-0.25 = 0.75 will be merged.
                                              deepSplit = 2, # medium sensitivity to cluster spliting
                                              minKMEtoStay = 0,
                                              pamRespectsDendro = FALSE, 
                                              minModuleSize = 30, 
                                              randomSeed = 90280, # there's some randomness associated with this calculation so we should set a seed
                                              nThreads = 6,
                                              verbose = 3)
    
    # Save coexpression network data   
    readr::write_rds(my.net.consensus, paste(dir.output, "TETranscripts_GEUVADIS_CONSENSUS_WGCNA_Network.RDS", sep = ""))
    
    # Extract network features for downstream use
    my.consMEs = my.net.consensus$multiMEs
    my.moduleColors = my.net.consensus$colors
    
    # Save the subset of network features
    save(my.consMEs, my.moduleColors, file = paste(dir.output, "TETranscripts_GEUVADIS_CONSENSUS_MEs_and_moduleColors.RData", sep = ''))
    

# Generate dendogram of the network modules
      
    # Extract gene module colors
    my.mergedColors.consensus = my.net.consensus[['colors']]

    # Plot the dendrogram and the module colors underneath
    pdf(paste(dir.output, "TETranscripts_GEUVADIS_Consensus_Module_Dendogram_Plot", ".pdf", sep=""))
    
        plotDendroAndColors(my.net.consensus$dendrograms[[1]], 
                            my.mergedColors.consensus[my.net.consensus$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, 
                            hang = 0.03,
                            addGuide = TRUE, 
                            guideHang = 0.05,
                            main = "Consensus gene dendrogram and module colors")
        
    dev.off()
    
    
# Make key linking genes to modules
      
    # Extract gene and module info  
    my.gene_module_key.consensus <- as.data.frame(tibble::enframe(my.net.consensus$colors, name = "gene", value = "module"))

    # Assign EnsemblID to rownames
    rownames(my.gene_module_key.consensus) <- my.gene_module_key.consensus$gene
    
    # Let's add the `ME` part so its more clear what these numbers are and it matches elsewhere
    my.gene_module_key.consensus$module <- paste("ME", my.gene_module_key.consensus$module, sep = "")

    # Add column specifying hub genes in each population
    my.gene_module_key.consensus$Hub_genes_EUR <- 'NA'
    my.gene_module_key.consensus$Hub_genes_YRI <- 'NA'

    # Identify hub genes
    my.hubs.EUR <- chooseTopHubInEachModule(datExpr = VST.expression.EUR[, my.gene_module_key.consensus$gene],
                                            colorh = my.gene_module_key.consensus$module,
                                            omitColors = c(""),
                                            power = 4,
                                            type = "signed")

    my.hubs.YRI <- chooseTopHubInEachModule(datExpr = VST.expression.YRI[, my.gene_module_key.consensus$gene],
                                            colorh = my.gene_module_key.consensus$module,
                                            omitColors = c(""),
                                            power = 4,
                                            type = "signed")

    # Fill in key table with hub gene designations (NOTE: only the top hub is shown, not multiple top hubs)
    my.gene_module_key.consensus[which(my.gene_module_key.consensus$gene %in% my.hubs.EUR), 'Hub_genes_EUR'] <- 'Hub_gene'
    my.gene_module_key.consensus[which(my.gene_module_key.consensus$gene %in% my.hubs.YRI), 'Hub_genes_YRI'] <- 'Hub_gene'
    
    # Add columns for alternative geneIDs (like gene symbols)
    my.gene_module_key.consensus <- Add_GeneSymbols_to_DF(input_df = my.gene_module_key.consensus, organism = 'hs')

    # Save module - gene key
    write.table(my.gene_module_key.consensus, file = paste(dir.output, "TETranscripts_GEUVADIS_Consensus_Module_vs_Gene_Key", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote=F)


# Make table/plot for the number of TEs per module
    
    # Keep only repeats
    key.repeats.consensus <- my.gene_module_key.consensus[!grepl('ENSG', my.gene_module_key.consensus$gene), ]

    # Add a column to hold counts
    key.repeats.consensus$counts <- 1

    # Aggregate the table by module, adding the counts
    repeats.per.module.consensus <- aggregate(x = key.repeats.consensus[, ncol(key.repeats.consensus)], by = list(key.repeats.consensus$module), FUN = 'sum')

    # Update the colnames
    colnames(repeats.per.module.consensus) <- c('Module', 'Counts')

    # Order from highest to lowest number of repeats
    repeats.per.module.consensus <- repeats.per.module.consensus[order(repeats.per.module.consensus$Counts, decreasing = TRUE), ]

    # Lock in factor level order
    repeats.per.module.consensus$Module <- factor(repeats.per.module.consensus$Module, levels = unique(repeats.per.module.consensus$Module))

    # Save results
    write.table(repeats.per.module.consensus, file = paste(dir.output, "TETranscripts_GEUVADIS_Consensus_Repeat_Counts_Per_Module", '.txt', sep =""), sep = "\t" , row.names = F, col.names = T, quote = F)

    # Generate bar plots
    pdf(paste(dir.output, "TETranscripts_GEUVADIS_Consensus_Repeat_Counts_Per_Module_Barplot", ".pdf", sep=""), width = 6, height = 6)
    par(mar=c(6,4,4,1)+.1)
    
        # Make a barplot
        barplot(repeats.per.module.consensus$Counts ~ repeats.per.module.consensus$Module,
                ylab = "# of Repeat Subfamilies",
                ylim = c(0, 1000),
                xlab = "",
                las = 2,
                cex.names = 1,
                main = 'Repeat Distribution per Cluster in Consensus Network'
                )
        
    dev.off()
    
    
    
    
    
# ASSESS EIGENGENE NETWORK PRESERVATION (EUR vs YRI)
    
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("European LCLs", "Yoruban LCLs")

# perform comparison and save the plot
pdf(file = paste(dir.output, "TETranscripts_GEUVADIS_CONSENSUS_Preservation_Across_EUR_and_YRI.pdf", sep = ''), width= 8, height = 10);
par(cex = 0.9)

    plotEigengeneNetworks(my.consMEs, 
                          plotAdjacency = FALSE,
                          setLabels, 
                          marDendro = c(0,2,2,1), 
                          marHeatmap = c(3,3,2,1),
                          zlimPreservation = c(0.5, 1), 
                          xLabelsAngle = 90,
                          excludeGrey = 'FALSE'
                          )
dev.off()
    
    



# MODULE EXPRESSION HEATMAPS ACROSS CONTROLS/CASES --------------------------------------------
         

# DEFINE INPUT PARAMETERS

# Define output directory
dir.output <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Trait_Correlations_Consensus/'

# Load network data 
load(file = '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Network_Construction_Consensus/TETranscripts_GEUVADIS_CONSENSUS_MEs_and_moduleColors.RData')

      # Extract module eigengene expression
      eigen.expr.EUR <- my.consMEs[[1]][["data"]]
      eigen.expr.YRI <- my.consMEs[[2]][["data"]]
      
      # transpose the data
      eigen.expr.EUR <- t(eigen.expr.EUR)
      eigen.expr.YRI <- t(eigen.expr.YRI)
  
# Load singleton counts
my.singletons <- read.csv("/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/1_Prepare_SV_Genotype_Data/Singleton_Frequency_and_Distribution_Analysis/Singleton_Counts_Per_Sample_L1_and_Alu.txt", header = T, stringsAsFactors = F, sep = '\t', row.names = 1)

    # Subset EUR and YRI
    singletons.EUR <- my.singletons[colnames(eigen.expr.EUR), ]
    singletons.YRI <- my.singletons[colnames(eigen.expr.YRI), ]
    
    # Order singleton samples by increasing singleton number
    singletons.EUR <- singletons.EUR[order(singletons.EUR$Combined_unique_insertions, decreasing = FALSE), ]
    singletons.YRI <- singletons.YRI[order(singletons.YRI$Combined_unique_insertions, decreasing = FALSE), ]
    
    # Add column with case/control label
    singletons.EUR$Label <- 'control'
    singletons.YRI$Label <- 'control'
    
    # Update the case labels
    singletons.EUR[which(singletons.EUR$Combined_unique_insertions > 0), 'Label'] <- 'case'
    singletons.YRI[which(singletons.YRI$Combined_unique_insertions > 0), 'Label'] <- 'case'
    
    # Update label column to a factor
    singletons.EUR$Label <- as.factor(singletons.EUR$Label)
    singletons.YRI$Label <- as.factor(singletons.YRI$Label)

# Reorder eigengene data by increasing singleton number
#eigen.expr.EUR <- eigen.expr.EUR[, rownames(singletons.EUR)]
#eigen.expr.YRI <- eigen.expr.YRI[, rownames(singletons.YRI)]

# Reorder eigengene data by MEroyalblue expression (correlated with case/control)

    # transpose the data
    eigen.expr.EUR <- as.data.frame(t(eigen.expr.EUR))
    eigen.expr.YRI <- as.data.frame(t(eigen.expr.YRI))
    
    # reorder
    eigen.expr.EUR <- eigen.expr.EUR[order(eigen.expr.EUR$MEroyalblue, decreasing = FALSE), ]
    eigen.expr.YRI <- eigen.expr.YRI[order(eigen.expr.YRI$MEroyalblue, decreasing = FALSE), ]
    
    # transpose back
    eigen.expr.EUR <- as.data.frame(t(eigen.expr.EUR))
    eigen.expr.YRI <- as.data.frame(t(eigen.expr.YRI))



    
    
# GENERATE HEATMAP (YRI)

# Define column label groups and group colors
annot_column <- data.frame(row.names = rownames(singletons.YRI), 
                           Plot_Group = singletons.YRI$Label)

annot_color <- list(Plot_Group = c(control = 'grey',
                                   case = 'black')
                    )
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)

# Save heatmap (DIFFERENTIAL EIGENGENES)
pdf(paste(dir.output, "Consensus_Network_EigenGene_Expression_YRI", ".pdf", sep=""), height = 8, width = 10)

pheatmap(eigen.expr.YRI, 
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         fontsize_row = 8,
         cluster_cols = FALSE, 
         show_colnames = FALSE,
         scale = 'row',
         border_color = NA,
         annotation_col = annot_column,
         annotation_colors = annot_color,
         color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
         breaks = breaksList
         #main = paste('bla \n', Number.Genes, ' Genes (FDR 5%)', sep = '')
         )
    
dev.off()





# GENERATE HEATMAP (EUR)

# Define column label groups and group colors
annot_column <- data.frame(row.names = rownames(singletons.EUR), 
                           Plot_Group = singletons.EUR$Label)

annot_color <- list(Plot_Group = c(control = 'grey',
                                   case = 'black')
                    )
  
# Define pheatmap scale
breaksList <- seq(-2, 2, by = 0.001)

# Save heatmap (DIFFERENTIAL EIGENGENES)
pdf(paste(dir.output, "Consensus_Network_EigenGene_Expression_EUR", ".pdf", sep=""), height = 8, width = 10)

pheatmap(eigen.expr.EUR, 
         cluster_rows = FALSE, 
         show_rownames = TRUE,
         fontsize_row = 8,
         cluster_cols = FALSE, 
         show_colnames = FALSE,
         scale = 'row',
         border_color = NA,
         annotation_col = annot_column,
         annotation_colors = annot_color,
         color = colorRampPalette(c("blue", "white", "red"))(length(breaksList)),
         breaks = breaksList
         #main = paste('bla \n', Number.Genes, ' Genes (FDR 5%)', sep = '')
         )
    
dev.off()






# Session Info ------------------------------------------------------------------------------------------------------------------------



# Save session info    
dir.session_info <- '/Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/2_RNASeq_Analyses/4_WGCNA/Session_Info/'
    
    
sink(file = paste(dir.session_info,"Session_Info_WGCNA_Network_Construction_Consensus.txt", sep =""))
sessionInfo()
sink()      
    


    
# Clean the environment
rm(list=ls())


