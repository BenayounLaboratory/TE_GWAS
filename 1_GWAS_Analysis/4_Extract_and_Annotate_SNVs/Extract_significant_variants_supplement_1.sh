# Go to directory with the significant SNVs
cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants

# Clump L1+Alu singleton GWAS SNVs using Plink
# Clump using the AFR panel, which has the most number of samples
# Make sure there are not duplicate RsIDs/SV_IDs, otherwise exclude them
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --bfile /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/Plink_AFR_Genotypes/Final_SNVs_SVs_AFR --clump Singleton_Meta-Analysis_Significant_Variants.txt --clump-p1 0.05 --clump-p2 1 --clump-r2 0.10 --clump-kb 250 --clump-snp-field SNP --clump-field 'P(R)' --set-missing-var-ids @:# --out Singleton_Meta-Analysis_Significant_Variants
