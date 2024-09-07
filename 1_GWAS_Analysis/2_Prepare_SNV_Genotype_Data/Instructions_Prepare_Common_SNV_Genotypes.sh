################### Prepare dbSNP build 155

# Go to directory where dbSNP data will be processed
cd SNV_Genotypes/dbSNP_Annotations

# Download the dbSNP VCF file (GCF_000001405.39.gz was used for this study), as well as the corresponding .tbi index file
# Link: https://ftp.ncbi.nih.gov/snp/latest_release/VCF/

# Download the corresponding genome assembly report (this will be used to map Refseq IDs in the dbSNP file to chromosome numbers)
# Link: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt

# Make a file with the chromosome refseq ids and chromosome number mappings
for k in *assembly_report.txt
  do
    out=$(echo $k | sed 's/.txt/.chrnames/')
    grep -e '^[^#]' $k | awk '{ print $7, $1 }' > $out
done

# Update the dbSNP file with chromosome numbers
bcftools annotate --rename-chrs GCF_000001405.39_GRCh38.p13_assembly_report.chrnames --threads 10 -Oz -o GRCh38.dbSNP155.vcf.gz GCF_000001405.39.gz


################### Annotate 1000Genomes SNVs with RsIDs

# Go up to the main genotype directory
cd ..

# Download GRCh38 1000Genomes Phase 3 SNV/Indel Data
# Link: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/supporting/related_samples/

# Index SNV files
for f in *.vcf.gz;do bcftools index $f --threads 10;done

# Apply dbSNP155 RsIDs
for f in *.vcf.gz;do bcftools annotate -a dbSNP_Annotations/GRCh38.dbSNP155.vcf.gz -c ID $f --threads 10 -Oz -o rsID.${f};done

# Delete or move the original SNV files elsewhere. dbSNP files can also be discarded or moved.

# Subset samples by super-population (FILTER 0)

    # Subset the samples. Manually move the output from each command to a separate folder. The VCF files with all samples can be deleted or moved afterwards.
    for f in *.vcf.gz;do bcftools view -S EUR_503_Both_Sexes -Oz -o EUR_503_Both.${f} $f --threads 10;done
    for f in *.vcf.gz;do bcftools view -S AFR_660_Both_Sexes -Oz -o AFR_660_Both.${f} $f --threads 10;done
    for f in *.vcf.gz;do bcftools view -S EAS_504_Both_Sexes -Oz -o EAS_504_Both.${f} $f --threads 10;done
    for f in *.vcf.gz;do bcftools view -S SAS_489_Both_Sexes -Oz -o SAS_489_Both.${f} $f --threads 10;done
    for f in *.vcf.gz;do bcftools view -S AMR_347_Both_Sexes -Oz -o AMR_347_Both.${f} $f --threads 10;done

    # Make a list of the vcf files for each superpopulation
    ls *vcf.gz > vcf_list

    # Manually reorganize each list so that it is ordered from chr1 to chr22 (autosomes)

    # Concatenate SNVs on each chromosome into one file
    bcftools concat --file-list vcf_list --output chr1-22.EUR_503_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive
    bcftools concat --file-list vcf_list --output chr1-22.AFR_660_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive
    bcftools concat --file-list vcf_list --output chr1-22.EAS_504_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive
    bcftools concat --file-list vcf_list --output chr1-22.SAS_489_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive
    bcftools concat --file-list vcf_list --output chr1-22.AMR_347_Both_rsID.phase3.genotypes.vcf.gz --threads 10 --naive

    # Update VCF headers for compatibility with vcftools
    sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.EUR_503_Both_rsID.phase3.genotypes.vcf > chr1-22.EUR_503_Both_rsID.4.2.phase3.genotypes.vcf
    sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.AFR_660_Both_rsID.phase3.genotypes.vcf > chr1-22.AFR_660_Both_rsID.4.2.phase3.genotypes.vcf
    sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.EAS_504_Both_rsID.phase3.genotypes.vcf > chr1-22.EAS_504_Both_rsID.4.2.phase3.genotypes.vcf
    sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.SAS_489_Both_rsID.phase3.genotypes.vcf > chr1-22.SAS_489_Both_rsID.4.2.phase3.genotypes.vcf
    sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' chr1-22.AMR_347_Both_rsID.phase3.genotypes.vcf > chr1-22.AMR_347_Both_rsID.4.2.phase3.genotypes.vcf

    # Move these concatenated files to the main 'SNV_Genotypes/Filtered_0' directory. Everything else can be moved or deleted.

    # Compress these files, but temporarily maintain a copy of uncompressed VCFs to run the next line of code


################### Filter analysis-grade SNVs with vcftools

# Filter common (MAF > 0.01) SNVs in each super-populations (FILTER #1)
vcftools --vcf chr1-22.EUR_503_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_chr1-22.EUR_503_BOTH_SEXES_chr1-22 --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --vcf chr1-22.AFR_660_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_chr1-22.AFR_660_BOTH_SEXES --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --vcf chr1-22.EAS_504_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_chr1-22.EAS_504_BOTH_SEXES --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --vcf chr1-22.SAS_489_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_chr1-22.SAS_489_BOTH_SEXES --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all
vcftools --vcf chr1-22.AMR_347_Both_rsID.4.2.phase3.genotypes.vcf --out Filtered_chr1-22.AMR_347_BOTH_SEXES --min-alleles 2 --max-alleles 2 --remove-indels --maf 0.01 --hwe .000001 --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all

    # Compress the files
    bgzip -@ 10 Filtered_chr1-22.EUR_503_BOTH_SEXES_chr1-22.recode.vcf
    bgzip -@ 10 Filtered_chr1-22.AFR_660_BOTH_SEXES.recode.vcf
    bgzip -@ 10 Filtered_chr1-22.EAS_504_BOTH_SEXES.recode.vcf
    bgzip -@ 10 Filtered_chr1-22.SAS_489_BOTH_SEXES.recode.vcf
    bgzip -@ 10 Filtered_chr1-22.AMR_347_BOTH_SEXES.recode.vcf

    # Index each VCF file
    for f in *.vcf.gz;do bcftools index $f --threads 10;done

    # For SNVs without an RsID, assign chr:pos for the ID
    bcftools annotate --set-id +'%CHROM\:%POS' Filtered_chr1-22.AFR_660_BOTH_SEXES.recode.vcf.gz > temporary.AFR.vcf
    bcftools annotate --set-id +'%CHROM\:%POS' Filtered_chr1-22.AMR_347_BOTH_SEXES.recode.vcf.gz > temporary.AMR.vcf
    bcftools annotate --set-id +'%CHROM\:%POS' Filtered_chr1-22.EAS_504_BOTH_SEXES.recode.vcf.gz > temporary.EAS.vcf
    bcftools annotate --set-id +'%CHROM\:%POS' Filtered_chr1-22.EUR_503_BOTH_SEXES_chr1-22.recode.vcf.gz > temporary.EUR.vcf
    bcftools annotate --set-id +'%CHROM\:%POS' Filtered_chr1-22.SAS_489_BOTH_SEXES.recode.vcf.gz > temporary.SAS.vcf

    # Make a list of common RsIDs in each super-population
    bcftools query -f '%ID\n' temporary.AFR.vcf > IDs.AFR.txt
    bcftools query -f '%ID\n' temporary.AMR.vcf > IDs.AMR.txt
    bcftools query -f '%ID\n' temporary.EAS.vcf > IDs.EAS.txt
    bcftools query -f '%ID\n' temporary.EUR.vcf > IDs.EUR.txt
    bcftools query -f '%ID\n' temporary.SAS.vcf > IDs.SAS.txt

    # Define common SNVs that are shared across all super-populations
    comm -12  <(sort IDs.AFR.txt) <(sort IDs.AMR.txt) | comm -12 - <(sort IDs.EAS.txt) | comm -12 - <(sort IDs.EUR.txt) | comm -12 - <(sort IDs.SAS.txt) > IDs.shared.txt


# Subset VCF files to include only common, shared SNVs (FILTER #2)
bcftools view -i '%ID=@IDs.shared.txt' temporary.AFR.vcf > Final_SNVs_AFR.recode.vcf
bcftools view -i '%ID=@IDs.shared.txt' temporary.AMR.vcf > Final_SNVs_AMR.recode.vcf
bcftools view -i '%ID=@IDs.shared.txt' temporary.EAS.vcf > Final_SNVs_EAS.recode.vcf
bcftools view -i '%ID=@IDs.shared.txt' temporary.EUR.vcf > Final_SNVs_EUR.recode.vcf
bcftools view -i '%ID=@IDs.shared.txt' temporary.SAS.vcf > Final_SNVs_SAS.recode.vcf

    # Compress VCF files
    bgzip -@ 10 Final_SNVs_AFR.recode.vcf
    bgzip -@ 10 Final_SNVs_AMR.recode.vcf
    bgzip -@ 10 Final_SNVs_EAS.recode.vcf
    bgzip -@ 10 Final_SNVs_EUR.recode.vcf
    bgzip -@ 10 Final_SNVs_SAS.recode.vcf

    # Index VCF files
    for f in *.vcf.gz;do bcftools index $f --threads 10;done

# Convert the final results from the AFR super-populations to PLINK BED format
/Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf Final_SNVs_AFR.recode.vcf.gz --set-missing-var-ids @:# --keep-allele-order --make-bed --out Plink_AFR_Genotypes/Final_SNVs_AFR
