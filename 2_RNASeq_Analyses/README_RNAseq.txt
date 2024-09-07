README -- RNAseq analyses
############################

1. Check for population structure in the SNV+SV genotypes among the GEUVADIS samples
    - Subset GEUVADIS samples, prune SNV+SV genotypes, and run PCA analysis (Instructions_Run_PCA_Analysis_GEUVADIS_samples_only.sh)
    - Plot SNV+SV genotype PCA results (Run_Population_Structure_Analysis_SNVs_and_SVs.R)

2. Prepare RNAseq data for WGCNA analysis
    - Filter lowly expressed genes, VST transform counts, and remove batch effects (Process_RNASeq.R, Process_RNASeq_functions.R)

3. Run DESeq2 analysis comparing (i) cases vs controls and (ii) genotypes across polymorphic structural variants (Run_DESeq.R, Run_DESeq_functions.R)

4. Construction and implementation of gene co-expression networks using WGCNA
    - Construct consensus gene co-expression networks using the EUR and AFR LCL batch-corrected expression data (Run_Consensus_WGCNA_Network_Construction.R, WGCNA_Functions.R)
    - Run network module-trait correlation analyses (Run_Consensus_WGCNA_Trait_Correlations.R, WGCNA_Functions.R)

5. Over-representation analysis (ORA) of consensus network modules
    - Prepare GMT gene sets for ORA (Prepare_Gene_Sets.R)
    - Run ORA for each module and plot, at most, the top 10 significant gene sets (Run_Overrepresentation_Analysis.R, Run_ORA_Functions.R)
