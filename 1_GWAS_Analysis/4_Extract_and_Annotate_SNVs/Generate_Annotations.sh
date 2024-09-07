# Change directory to raw results folder
cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations

# Intersect all GWAS SNVs (ONLY) and all SVs (ONLY) with STRUCTURAL VARIANT HOTSPOT (Lin and Gokcumen) bed file. Note, chromosomes should be in 'chr' format
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Lin_and_Gokcumen_2019_Structural_Variation_Hotspots/hg38_SV_Hotspots_Autosomal_Mapped_with_UCSC_Lift_Genome_Annotations_on_May_20_2024.bed -wa -wb > SNV_Variants_Intersecting_SV-Hotspots_Lin-2019.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Lin_and_Gokcumen_2019_Structural_Variation_Hotspots/hg38_SV_Hotspots_Autosomal_Mapped_with_UCSC_Lift_Genome_Annotations_on_May_20_2024.bed -wa -wb > SV_Variants_Intersecting_SV-Hotspots_Lin-2019.txt

# Intersect all GWAS SNVs (ONLY) and all SVs (ONLY) with SEGMENTAL DUPLICATION bed file (Bailey 2001). Note, chromosomes should be in 'chr' format
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bailey_2001_Segmental_Duplications/UCSC_hg38_Segmental_Duplicates_Updated_2014_10_14.bed -wa -wb > SNV_Variants_Intersecting_Segmental_Duplications_Bailey-2001.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/Bailey_2001_Segmental_Duplications/UCSC_hg38_Segmental_Duplicates_Updated_2014_10_14.bed -wa -wb > SV_Variants_Intersecting_Segmental_Duplications_Bailey-2001.txt

# Intersect all GWAS SNVs (ONLY) and all SVs (ONLY) with Repeatmasker bed file. Note, chromosomes should be in 'chr' format
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/UCSC_Repeatmasker_Annotations/Repeatmasker_chr_May_16_2022.bed -wa -wb > SNV_Variants_Intersecting_Repeatmasker.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/UCSC_Repeatmasker_Annotations/Repeatmasker_chr_May_16_2022.bed -wa -wb > SV_Variants_Intersecting_Repeatmasker.txt

# Intersect all GWAS SNVs (ONLY) and all SVs (ONLY) with ENCODE Blacklist v2 bed file. Note, chromosomes should be in 'chr' format
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Blacklists/hg38-blacklist.v2.bed -wa -wb > SNV_Variants_Intersecting_ENCODE_Blacklist.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Blacklists/hg38-blacklist.v2.bed -wa -wb > SV_Variants_Intersecting_ENCODE_Blacklist.txt

# Intersect all GWAS SNVs with L1Base2 bed file for full length, intact-for-both-ORFs L1s
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/L1Base_V2/NCBI38_chr_Human_Full_Length_Intact_L1_06_03_2022_hsflil1_8438.bed -wa -wb > SNV_Variants_Intersecting_L1Base2_Full_Length_L1s.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/L1Base_V2/NCBI38_chr_Human_Full_Length_Intact_L1_06_03_2022_hsflil1_8438.bed -wa -wb > SV_Variants_Intersecting_L1Base2_Full_Length_L1s.txt

# Intersect all GWAS SNVs with L1Base2 bed file for disrupted ORF1 but intact ORF2 L1s
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/L1Base_V2/NCBI38_chr_Human_ORF2_Intact_L1_06_03_2022_hsorf2l1_8438.bed -wa -wb > SNV_Variants_Intersecting_L1Base2_Intact_ORF2_L1.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/L1Base_V2/NCBI38_chr_Human_ORF2_Intact_L1_06_03_2022_hsorf2l1_8438.bed -wa -wb > SV_Variants_Intersecting_L1Base2_Intact_ORF2_L1.txt

# Intersect all GWAS SNVs with L1Base2 bed file for full length, non-intact L1s
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/L1Base_V2/NCBI38_chr_Human_Full_Length_NonIntact_L1_06_03_2022_hsflnil1_8438_rm.bed -wa -wb > SNV_Variants_Intersecting_L1Base2_NonIntact_Full_Length_L1s.txt
bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/L1Base_V2/NCBI38_chr_Human_Full_Length_NonIntact_L1_06_03_2022_hsflnil1_8438_rm.bed -wa -wb > SV_Variants_Intersecting_L1Base2_NonIntact_Full_Length_L1s.txt

# Split all GWAS SNVs into files with 1 million SNVs and submit each to the GREAT v4.0.4 regulatory annotation platform, using these settings: GRCh38 (UCSC hg38, Dec. 2013), whole genome background, basal plus extension, 5 kb upstream and 1 kb downstream, plus distal up to 1000 kb
# No need to split SVs since there are only several thousand; submit bed file to the GREAT v4.0.4 regulatory annotation platform, using these settings: GRCh38 (UCSC hg38, Dec. 2013), whole genome background, basal plus extension, 5 kb upstream and 1 kb downstream, plus distal up to 1000 kb

# Intersect all GWAS SNVs with the ENCODE Registry 4 of cCREs
cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/ENCODE_cCRE_Annotations

    # FOR SNVs
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.bed -wa -wb > SNV_Variants_Intersecting_ENCODE_cCREs.txt
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.CTCF-only.bed -wa -wb > SNV_Variants_Intersecting_ENCODE_CTCF_Bound.txt
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.ELS.bed -wa -wb > SNV_Variants_Intersecting_ENCODE_Candidate_Enhancers.txt
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SNVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.PLS.bed -wa -wb > SNV_Variants_Intersecting_ENCODE_Candidate_Promoters.txt

    # FOR SVs
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.bed -wa -wb > SV_Variants_Intersecting_ENCODE_cCREs.txt
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.CTCF-only.bed -wa -wb > SV_Variants_Intersecting_ENCODE_CTCF_Bound.txt
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.ELS.bed -wa -wb > SV_Variants_Intersecting_ENCODE_Candidate_Enhancers.txt
    bedtools intersect -a /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Resource_Files_SNVs_and_SVs/RESOURCE_All_Background_Positions_SVs_only.bed -b /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/0_Sample_Metadata_and_External_Resources/ENCODE_Registry_v4_hg38_cCREs/GRCh38-cCREs.PLS.bed -wa -wb > SV_Variants_Intersecting_ENCODE_Candidate_Promoters.txt

# ANNOTATE WITH SNPEFF/SNPSIFT

    # PREPARE VCF FILE

        # Change to the target directory
        cd /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Raw_Annotations/SNPEff

        # Keep only interesting variants in VCF format
        bcftools view -T /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/4_Extract_and_Annotate_SNVs/Unannotated_Significant_Variants/Singleton_Meta-Analysis_Significant_Variants_Unique_Variant_Positions.txt -Oz -o Combined_Singleton_Meta-Analysis_Significant_Variants.vcf.gz /Users/juanb/Desktop/2024_TE_GWAS_Juan/Code/1_GWAS_Analysis/2_Prepare_SNV_Genotype_Data/Combined_SNV_SV_Genotypes/AFR.SNVs.SVs.vcf.gz --threads 8

    # RUN SNPEFF/SNPSIFT

        # base SNPEff annotations
        java -Xmx32g -jar /Users/juanb/Documents/Bioinformatic_Tools/snpEff_5.2c/snpEff.jar -v -motif GRCh38.99 Combined_Singleton_Meta-Analysis_Significant_Variants.vcf.gz > Combined_Singleton_Meta-Analysis_Significant_Variants.SNPEff.vcf

        # NOTE: PhastCons files can be downloaded here
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr1.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr2.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr3.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr4.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr5.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr6.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr7.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr8.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr9.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr10.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr11.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr12.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr13.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr14.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr15.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr16.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr17.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr18.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr19.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr20.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr21.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr22.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chrM.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chrX.phastCons100way.wigFix.gz
        # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chrY.phastCons100way.wigFix.gz

        # phastCons annotation
        java -Xmx32g -jar /Users/juanb/Documents/Bioinformatic_Tools/snpEff_5.2c/SnpSift.jar phastCons /Users/juanb/Documents/Bioinformatic_Tools/snpEff_5.2c/db/phastCons_hg38 Combined_Singleton_Meta-Analysis_Significant_Variants.SNPEff.vcf > Combined_Singleton_Meta-Analysis_Significant_Variants.SNPEff.phastCons.vcf

        # Extract annotations in an easy-to-read format
        cat Combined_Singleton_Meta-Analysis_Significant_Variants.SNPEff.phastCons.vcf | /Users/juanb/Documents/Bioinformatic_Tools/snpEff_5.2c/scripts/vcfEffOnePerLine.pl | java -jar /Users/juanb/Documents/Bioinformatic_Tools/snpEff_5.2c/SnpSift.jar extractFields - CHROM POS ID REF ALT "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATUREID" "ANN[*].HGVS_P" "PhastCons" > READABLE_ANNOTATIONS.csv
