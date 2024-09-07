# Download Phase 3 structural variant calls and the corresponding .tbi index file
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/GRCh38_positions/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz

# Extract autosomal global singleton SVs, using VCFtools, and compress the output using bgzip
vcftools --gzvcf ALL.wgs.mergedSV.v8.20130502.svs.genotypes.GRCh38.vcf.gz --out Single_ALT_alleles.wgs.mergedSV.v8.20130502.svs.genotypes --max-missing 1 --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --recode --recode-INFO-all --non-ref-ac 1 --max-non-ref-ac 1
bgzip -@ 10 Single_ALT_alleles.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf

# Extract L1 and Alu insertion global singletons, using bcftools
bcftools view -i 'SVTYPE ="LINE1"' Single_ALT_alleles.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf --threads 10 -Oz -o ALL_unique_LINE1.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
bcftools view -i 'SVTYPE ="ALU"' Single_ALT_alleles.wgs.mergedSV.v8.20130502.svs.genotypes.recode.vcf --threads 10 -Oz -o ALL_unique_ALU.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Extract L1 and Alu insertion global singletons that contain target site duplications (TSDs) (NOT USED IN THE FINAL STUDY)
bcftools view -e 'TSD ="null"' ALL_unique_LINE1.wgs.mergedSV.v8.20130502.svs.genotypes.vcf --threads 10 -Oz -o With_TSD_ALL_unique_LINE1.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz
bcftools view -e 'TSD ="null"' ALL_unique_ALU.wgs.mergedSV.v8.20130502.svs.genotypes.vcf --threads 10 -Oz -o With_TSD_ALL_unique_ALU.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

# Convert all singleton VCFs to PLINK BED format

    # With and without TSDs
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf ALL_unique_LINE1.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out ALL_unique_LINE1.recode.plink
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf ALL_unique_ALU.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out ALL_unique_ALU.recode.plink

    # Only with TSDs
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf With_TSD_ALL_unique_LINE1.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out With_TSD_ALL_unique_LINE1.recode.plink
    /Users/juanb/Documents/Bioinformatic_Tools/plink_mac_20200428/plink --vcf With_TSD_ALL_unique_ALU.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz --make-bed --out With_TSD_ALL_unique_ALU.recode.plink
