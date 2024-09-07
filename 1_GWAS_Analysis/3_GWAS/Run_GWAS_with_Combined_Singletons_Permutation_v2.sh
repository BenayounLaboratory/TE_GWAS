# Loop over 20 permutations

for permutation_number in {1..20}; do

  # Make directory to hold results for the ith permutation
  mkdir /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number

  # Run the ith permutation for AFR samples
  cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes
  /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_AFR --keep-allele-order --logistic --sex --covar Covariates_AFR_660.txt --pheno /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/AFR/Phenotype_Permutation_$permutation_number.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number/GWAS_Combined_Singletons_AFR_$permutation_number --hide-covar --threads 10

  # Run the ith permutation for AMR samples
  cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_Genotypes
  /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_AMR --keep-allele-order --logistic --sex --covar Covariates_AMR_347.txt --pheno /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/AMR/Phenotype_Permutation_$permutation_number.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number/GWAS_Combined_Singletons_AMR_$permutation_number --hide-covar --threads 10

  # Run the ith permutation for EAS samples
  cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_Genotypes
  /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_EAS --keep-allele-order --logistic --sex --covar Covariates_EAS_504.txt --pheno /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/EAS/Phenotype_Permutation_$permutation_number.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number/GWAS_Combined_Singletons_EAS_$permutation_number --hide-covar --threads 10

  # Run the ith permutation for EUR samples
  cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes
  /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_EUR --keep-allele-order --logistic --sex --covar Covariates_EUR_503.txt --pheno /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/EUR/Phenotype_Permutation_$permutation_number.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number/GWAS_Combined_Singletons_EUR_$permutation_number --hide-covar --threads 10

  # Run the ith permutation for SAS samples
  cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_Genotypes
  /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_SAS --keep-allele-order --logistic --sex --covar Covariates_SAS_489.txt --pheno /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permuted_Phenotypes/SAS/Phenotype_Permutation_$permutation_number.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number/GWAS_Combined_Singletons_SAS_$permutation_number --hide-covar --threads 10

  # Run the ith meta-analysis
  cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_Permutations/Permutation_$permutation_number
  /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --meta-analysis GWAS_Combined_Singletons_AFR_$permutation_number.assoc.logistic GWAS_Combined_Singletons_AMR_$permutation_number.assoc.logistic GWAS_Combined_Singletons_EAS_$permutation_number.assoc.logistic GWAS_Combined_Singletons_EUR_$permutation_number.assoc.logistic GWAS_Combined_Singletons_SAS_$permutation_number.assoc.logistic

done
