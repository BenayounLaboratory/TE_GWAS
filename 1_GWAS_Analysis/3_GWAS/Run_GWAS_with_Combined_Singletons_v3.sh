# GWAS with ALL global Alu+L1 singletons

    # AFR
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_AFR --keep-allele-order --logistic --sex --covar Covariates_AFR_660.txt --pheno Phenotype_Combined_L1_Alu_Singleton_AFR_660.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AFR_GWAS --hide-covar --threads 12

    # AMR
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_AMR --keep-allele-order --logistic --sex --covar Covariates_AMR_347.txt --pheno Phenotype_Combined_L1_Alu_Singleton_AMR_347.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/AMR_GWAS --hide-covar --threads 12

    # EAS
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_EAS --keep-allele-order --logistic --sex --covar Covariates_EAS_504.txt --pheno Phenotype_Combined_L1_Alu_Singleton_EAS_504.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EAS_GWAS --hide-covar --threads 12

    # EUR
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_EUR --keep-allele-order --logistic --sex --covar Covariates_EUR_503.txt --pheno Phenotype_Combined_L1_Alu_Singleton_EUR_503.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/EUR_GWAS --hide-covar --threads 12

    # SAS
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_SAS --keep-allele-order --logistic --sex --covar Covariates_SAS_489.txt --pheno Phenotype_Combined_L1_Alu_Singleton_SAS_489.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results/SAS_GWAS --hide-covar --threads 12

    # META Analysis
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --meta-analysis AFR_GWAS.assoc.logistic AMR_GWAS.assoc.logistic EAS_GWAS.assoc.logistic EUR_GWAS.assoc.logistic SAS_GWAS.assoc.logistic

# Subset GWAS with global Alu+L1 singletons CONTAINING A TARGET SITE DUPLICATION (THIS IS EXPLORATORY)

    # AFR
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_AFR --keep-allele-order --logistic --sex --covar Covariates_AFR_660.txt --pheno Phenotype_TSD_Combined_L1_Alu_Singleton_AFR_660.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_TSDs/AFR_GWAS --hide-covar --threads 12

    # AMR
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AMR_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_AMR --keep-allele-order --logistic --sex --covar Covariates_AMR_347.txt --pheno Phenotype_TSD_Combined_L1_Alu_Singleton_AMR_347.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_TSDs/AMR_GWAS --hide-covar --threads 12

    # EAS
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EAS_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_EAS --keep-allele-order --logistic --sex --covar Covariates_EAS_504.txt --pheno Phenotype_TSD_Combined_L1_Alu_Singleton_EAS_504.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_TSDs/EAS_GWAS --hide-covar --threads 12

    # EUR
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_EUR_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_EUR --keep-allele-order --logistic --sex --covar Covariates_EUR_503.txt --pheno Phenotype_TSD_Combined_L1_Alu_Singleton_EUR_503.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_TSDs/EUR_GWAS --hide-covar --threads 12

    # SAS
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_SAS_Genotypes
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile Final_SNVs_SVs_SAS --keep-allele-order --logistic --sex --covar Covariates_SAS_489.txt --pheno Phenotype_TSD_Combined_L1_Alu_Singleton_SAS_489.txt --ci 0.95 --out /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_TSDs/SAS_GWAS --hide-covar --threads 12

    # META Analysis
    cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/3_GWAS/GWAS_Raw_Results_TSDs
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --meta-analysis AFR_GWAS.assoc.logistic AMR_GWAS.assoc.logistic EAS_GWAS.assoc.logistic EUR_GWAS.assoc.logistic SAS_GWAS.assoc.logistic
