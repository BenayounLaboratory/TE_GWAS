################### Prepare the final SNV+SV genotype files

# For each super-population, make a file that holds the file names for the SNV and SV genotype files

# Combine the SNV and SV genotype files
bcftools concat --file-list file_list_AFR.txt --output AFR.SNVs.SVs.vcf --threads 10 -a
bcftools concat --file-list file_list_AMR.txt --output AMR.SNVs.SVs.vcf --threads 10 -a
bcftools concat --file-list file_list_EAS.txt --output EAS.SNVs.SVs.vcf --threads 10 -a
bcftools concat --file-list file_list_EUR.txt --output EUR.SNVs.SVs.vcf --threads 10 -a
bcftools concat --file-list file_list_SAS.txt --output SAS.SNVs.SVs.vcf --threads 10 -a

    # Compress the files
    bgzip -@ 10 AFR.SNVs.SVs.vcf
    bgzip -@ 10 AMR.SNVs.SVs.vcf
    bgzip -@ 10 EAS.SNVs.SVs.vcf
    bgzip -@ 10 EUR.SNVs.SVs.vcf
    bgzip -@ 10 SAS.SNVs.SVs.vcf


# SNV+SV genotype PCA Analysis with Plink

    # Prune genotypes in each super-population
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf AFR.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_AFR_PCA/Final_SNV_SV_AFR_pruned
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf AMR.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_AMR_PCA/Final_SNV_SV_AMR_pruned
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf EAS.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_EAS_PCA/Final_SNV_SV_EAS_pruned
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf EUR.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_EUR_PCA/Final_SNV_SV_EUR_pruned
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf SAS.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out Plink_SAS_PCA/Final_SNV_SV_SAS_pruned

    # Use the pruned data to run PCA analysis in each super-population
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf AFR.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_AFR_PCA/Final_SNV_SV_AFR_pruned.prune.in --pca var-wts --out Plink_AFR_PCA/Final_SNV_SV_AFR_PCA
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf AMR.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_AMR_PCA/Final_SNV_SV_AMR_pruned.prune.in --pca var-wts --out Plink_AMR_PCA/Final_SNV_SV_AMR_PCA
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf EAS.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_EAS_PCA/Final_SNV_SV_EAS_pruned.prune.in --pca var-wts --out Plink_EAS_PCA/Final_SNV_SV_EAS_PCA
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf EUR.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_EUR_PCA/Final_SNV_SV_EUR_pruned.prune.in --pca var-wts --out Plink_EUR_PCA/Final_SNV_SV_EUR_PCA
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf SAS.SNVs.SVs.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract Plink_SAS_PCA/Final_SNV_SV_SAS_pruned.prune.in --pca var-wts --out Plink_SAS_PCA/Final_SNV_SV_SAS_PCA

# Make the necessary covariate and phenotype files in R (i.e. run Run_Population_Structure_Analysis_SNVs_and_SVs.R and SNV_Genotype_Preparation_Supplement_1.R)

# Generate Plink BED files for the final SNV+SV genotypes
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf AFR.SNVs.SVs.vcf.gz --set-missing-var-ids @:# --keep-allele-order --keep Plink_AFR_Genotypes/Keep_AFR_660.txt --update-sex Plink_AFR_Genotypes/Sex_AFR_660.txt --make-bed --out Plink_AFR_Genotypes/Final_SNVs_SVs_AFR
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf AMR.SNVs.SVs.vcf.gz --set-missing-var-ids @:# --keep-allele-order --keep Plink_AMR_Genotypes/Keep_AMR_347.txt --update-sex Plink_AMR_Genotypes/Sex_AMR_347.txt --make-bed --out Plink_AMR_Genotypes/Final_SNVs_SVs_AMR
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf EAS.SNVs.SVs.vcf.gz --set-missing-var-ids @:# --keep-allele-order --keep Plink_EAS_Genotypes/Keep_EAS_504.txt --update-sex Plink_EAS_Genotypes/Sex_EAS_504.txt --make-bed --out Plink_EAS_Genotypes/Final_SNVs_SVs_EAS
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf EUR.SNVs.SVs.vcf.gz --set-missing-var-ids @:# --keep-allele-order --keep Plink_EUR_Genotypes/Keep_EUR_503.txt --update-sex Plink_EUR_Genotypes/Sex_EUR_503.txt --make-bed --out Plink_EUR_Genotypes/Final_SNVs_SVs_EUR
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf SAS.SNVs.SVs.vcf.gz --set-missing-var-ids @:# --keep-allele-order --keep Plink_SAS_Genotypes/Keep_SAS_489.txt --update-sex Plink_SAS_Genotypes/Sex_SAS_489.txt --make-bed --out Plink_SAS_Genotypes/Final_SNVs_SVs_SAS

# Make a list of duplicate IDs (this will be useful for clumping SNVs downstream). This only needs to be run once since all samples use the same SNV/SV list
cut -f 2 Plink_AFR_Genotypes/Final_SNVs_SVs_AFR.bim | sort | uniq -d > SNV_SV_duplicate_IDs.txt
