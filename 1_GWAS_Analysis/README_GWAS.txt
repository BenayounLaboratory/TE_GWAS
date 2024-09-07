README -- GWAS analyses
############################

1. Prepare structural variant (SV) genotype data for GWAS analysis
    - Extract L1 and Alu insertion global singletons, which will be used as the phenotype in the GWAS (Instructions_Extract_Global_Singletons.sh)
    - Extract common and shared (across all 5 super-populations) SVs to include in our GWAS scan (Instructions_Prepare_Common_SV_Genotypes.sh)
    - Process the global singleton data and generate various frequency/distribution plots (Singleton_Quantification_Supplement.R)

2. Prepare single nucleotide variant (SNV) genotype data for GWAS
    - Annotate variants with dbSNP RsIDs and extract common and shared SNVs to include in our GWAS scan (Instructions_Prepare_Common_SNV_Genotypes.sh)
    - Generate the final SNV+SV genotype matrices to use in the GWAS, and prune genotypes for PCA analysis (Instructions_Prepare_Final_Combined_SNVs_and_SV_Genotypes.sh)
    - Plot SNV+SV genotype PCA results (Run_Population_Structure_Analysis_SNVs_and_SVs.R)
    - Prepare covariate and phenotype data (i.e. presence or absense of L1 or Alu singletons) for PLINK (SNV_Genotype_Preparation_Supplement_1.R)
    - Generate SNV and SV ID/allele/position/etc. mapping files (SNV_Genotype_Preparation_Supplement_2.R)

3. Run all GWAS analyses
    - Take all sample-phenotype mappings, scramble the samples, and generate files with scrambled sample-phenotype mappings (Generate_Permuted_Samples.R, Generate_Permuted_Samples_Functions.R)
    - Run GWAS in each super-population and carry out the meta-analysis (Run_GWAS_with_Combined_Singletons_v3.sh)
    - Carry out 20 GWAS analyses in each super-population using the scrambled phenotypes, and use these results to generate 20 scrambled meta-analyses (Run_GWAS_with_Combined_Singletons_Permutation_v2.sh)

4. Extract and annotate significant SNVs and SVs
    - Calculate empirical FDR thresholds using the permutation results and extract * unannotated * significant variants (Extract_significant_variants.R, Extract_significant_variants_functions.R)
    - Clump significant variants using the AFR sample panel (Extract_significant_variants_supplement_1.sh)
    - Generate * raw * variant-region annotations (all variants) and significant variant-SnpEff annotations (Generate_Annotations.sh)
        -NOTE: All variants were also submitted to GREAT to generate raw variant-gene annotations.
    - Apply annotations to all variants and subset * annotated * significant variants. Add any remaining significant variant-specific annotations (Annotate_SNVs.R)

5. Determine whether variants are enriched for particular annotations
    - Determine whether significant SNVs are enriched in particular regions or for particular gene sets (Calculate_Enrichments_SNVs.R, Calculate_Enrichments_SNVs_functions.R)
    - Run over-representation analysis on SNV and SV-associated genes (Run_Overrepresentation_Analysis.R, Run_ORA_Functions.R)
    - Generate a heatmap for significant SNV-associated gene expression, in the GTEx Analysis v8  (GTEx_Analysis.R)
    - Run online enrichment tools
        - TEENA: Run on Aug 8 2024 with these settings: hg38 genome assembly, (yes) use mid-point for overlapping, (yes) exclude genomic gaps, and (no) exclude promoter regions
        - GREAT v4.0.4: Run on Aug 8 2024 with these settings for SNVs: GRCh38 (UCSC hg38, Dec. 2013), Whole Genome Background Regions, Basal Plus Extension rule with (local) 5 kb upstream, 1 kb downstream plus (distal) up to 1000 kb
                      - For SVs: Use all SVs that made it to the GWAS as the background, not the whole genome, since SVs have much less coverage of the entire genome.

6. Generate various GWAS plots
    - Generate manhattan plots for the GWAS results (Generate_Manhattan_Plots.R, Generate_Manhattan_Plots_Functions.R)
    - For interesting variants, generate beeswarm plots comparing genotypes vs fractions of samples with a global L1/Alu insertion singleton (Plot_Singletons_vs_Genotypes.R, Plot_Singletons_vs_Genotypes_Functions.R)
