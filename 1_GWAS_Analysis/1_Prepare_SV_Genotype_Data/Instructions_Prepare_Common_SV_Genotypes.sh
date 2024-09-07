# Filter common (MAF > 0.01) SVs in each super-populations (FILTER #1)
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz --out SVs_FILTERED_AFR_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_AFR_660 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz --out SVs_FILTERED_AMR_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_AMR_347 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz --out SVs_FILTERED_EAS_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_EAS_504 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz --out SVs_FILTERED_EUR_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_EUR_503 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz --out SVs_FILTERED_SAS_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --maf 0.01 --hwe .000001 --max-missing 1 --keep /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_SAS_489 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all

# Organize samples in a specific order, and use that same order to organize SNVs later on (necessary for concatenating SVs and SNVs)
bcftools view -S /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_AFR_660 -o reordered_SVs_FILTERED_AFR_BOTH_SEXES_chr1-22.recode.vcf SVs_FILTERED_AFR_BOTH_SEXES_chr1-22.recode.vcf
bcftools view -S /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_AMR_347 -o reordered_SVs_FILTERED_AMR_BOTH_SEXES_chr1-22.recode.vcf SVs_FILTERED_AMR_BOTH_SEXES_chr1-22.recode.vcf
bcftools view -S /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_EAS_504 -o reordered_SVs_FILTERED_EAS_BOTH_SEXES_chr1-22.recode.vcf SVs_FILTERED_EAS_BOTH_SEXES_chr1-22.recode.vcf
bcftools view -S /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_EUR_503 -o reordered_SVs_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf SVs_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf
bcftools view -S /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/1000G_Samples/Both_Sexes_SAS_489 -o reordered_SVs_FILTERED_SAS_BOTH_SEXES_chr1-22.recode.vcf SVs_FILTERED_SAS_BOTH_SEXES_chr1-22.recode.vcf

    # Compress files
    bgzip -@ 10 reordered_SVs_FILTERED_AFR_BOTH_SEXES_chr1-22.recode.vcf
    bgzip -@ 10 reordered_SVs_FILTERED_AMR_BOTH_SEXES_chr1-22.recode.vcf
    bgzip -@ 10 reordered_SVs_FILTERED_EAS_BOTH_SEXES_chr1-22.recode.vcf
    bgzip -@ 10 reordered_SVs_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf
    bgzip -@ 10 reordered_SVs_FILTERED_SAS_BOTH_SEXES_chr1-22.recode.vcf

    # Index each file
    for f in *.vcf.gz;do bcftools index $f --threads 10;done

# Make a list of common SVs in each super-population
bcftools query -f '%ID\n' reordered_SVs_FILTERED_AFR_BOTH_SEXES_chr1-22.recode.vcf.gz > IDs.AFR.txt
bcftools query -f '%ID\n' reordered_SVs_FILTERED_AMR_BOTH_SEXES_chr1-22.recode.vcf.gz > IDs.AMR.txt
bcftools query -f '%ID\n' reordered_SVs_FILTERED_EAS_BOTH_SEXES_chr1-22.recode.vcf.gz > IDs.EAS.txt
bcftools query -f '%ID\n' reordered_SVs_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz > IDs.EUR.txt
bcftools query -f '%ID\n' reordered_SVs_FILTERED_SAS_BOTH_SEXES_chr1-22.recode.vcf.gz > IDs.SAS.txt

# Define common SVs that are shared across all super-populations
comm -12  <(sort IDs.AFR.txt) <(sort IDs.AMR.txt) | comm -12 - <(sort IDs.EAS.txt) | comm -12 - <(sort IDs.EUR.txt) | comm -12 - <(sort IDs.SAS.txt) > IDs.shared.txt

# Subset VCF files to include only common, shared SVs (FILTER #2)
bcftools view -i ID=@IDs.shared.txt reordered_SVs_FILTERED_AFR_BOTH_SEXES_chr1-22.recode.vcf.gz > Final_SVs_AFR.recode.vcf
bcftools view -i ID=@IDs.shared.txt reordered_SVs_FILTERED_AMR_BOTH_SEXES_chr1-22.recode.vcf.gz > Final_SVs_AMR.recode.vcf
bcftools view -i ID=@IDs.shared.txt reordered_SVs_FILTERED_EAS_BOTH_SEXES_chr1-22.recode.vcf.gz > Final_SVs_EAS.recode.vcf
bcftools view -i ID=@IDs.shared.txt reordered_SVs_FILTERED_EUR_BOTH_SEXES_chr1-22.recode.vcf.gz > Final_SVs_EUR.recode.vcf
bcftools view -i ID=@IDs.shared.txt reordered_SVs_FILTERED_SAS_BOTH_SEXES_chr1-22.recode.vcf.gz > Final_SVs_SAS.recode.vcf

    # Compress VCF files
    bgzip -@ 10 Final_SVs_AFR.recode.vcf
    bgzip -@ 10 Final_SVs_AMR.recode.vcf
    bgzip -@ 10 Final_SVs_EAS.recode.vcf
    bgzip -@ 10 Final_SVs_EUR.recode.vcf
    bgzip -@ 10 Final_SVs_SAS.recode.vcf

    # Index VCF files
    for f in *.vcf.gz;do bcftools index $f --threads 10;done

# Convert the final results from the AFR super-populations to PLINK BED format
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Final_SVs_AFR.recode.vcf.gz --set-missing-var-ids @:# --keep-allele-order --make-bed --out Plink_AFR_Genotypes/Final_SVs_AFR
