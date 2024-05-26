#!/bin/bash

# Over script

mkdir -p ~/variant_calling/{aligned_reads,data,reads,scripts,results,tools,referentie_genoom}

# Download data

wget -P ~/variant_calling/reads ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz

wget -P ~/variant_calling/reads ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


echo "Run Prep files..."

################################ Prep files (to Be Generated Only Once) ###################################################

# Download reference files

wget -P ~/variant_calling/referentie_genoom/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/variant_calling/referentie_genoom/hg38.fa.gz

# Conda ENV maken

conda creat -n variant
conda activate variant

# Samtools installeren
#conda install -c bioconda samtools

# index ref- .fai file before running haplotype caller

#samtools faidx ~/variant_calling/referentie_genoom/hg38.fa

# GATK downloden en installeren
wget -P ~/variant_calling/tools https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
cd ~/variant_calling/tools
unzip gatk-4.5.0.0.zip
cd

# java installeren met conda
conda install -c conda-forge openjdk=17

# ref dict -.dict file before running haplotype caller
java -jar ~/variant_calling/tools/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CreateSequenceDictionary R=~/variant_calling/referentie_genoom/hg38.fa O=~/variant_calling/referentie_genoom/hg38.dict

# Download know sites files for BQSR from GATK resource bundle

wget -P ~/variant_calling/referentie_genoom/ https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ~/variant_calling/referentie_genoom/ https://storage.googleapis.com/genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

###################### Variant Calling Steps ###################################

# directories 

ref="$HOME/variant_calling/referentie_genoom/hg38.fa"
know_sites="$HOME/variant_calling/referentie_genoom/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="$HOME/variant_calling/aligned_reads"
reads="$HOME/variant_calling/reads"
results="$HOME/variant_calling/results"
tools="$HOME/variant_calling/tools"
data="$HOME/variant_calling/data"


#------------------------------
# Stap 1: QC met fastqc runnen
#------------------------------
echo "Stap1: QC met fastqc runnen"

 fastqc installeren met conda
conda install -c bioconda fastqc


fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No tremmimg required, quality looks okay.

#----------------------------------------
# Stap 2: Map to reference using BWA-MEM
#----------------------------------------

echo "Stap 2: Map to reference using BWA-MEM"

# BWA installeren met conda
conda install -c bioconda bwa

# BWA index reference
bwa index ${ref}

# BWA aligment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam


#------------------------------------------
# Stap 3: Mark Duplicates and Sort - GATK4
#------------------------------------------

echo "Stap 3: Mark Duplicates and Sort - GATK4"
java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O  ${aligned_reads}/SRR062634_sorted_dedup_reads.bam


#------------------------------------
# Stap 4: Base quality recalibration
#------------------------------------

echo "Stap 4: Base quality recalibration"

# 1. build the model
java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${know_sites} -O ${data}/recal_data.table

# 2. Apply the model to adjust the base quality scores
java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam

# Stap 5: Collect Alignment & insert Size Metrics

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf


#----------------------------------------------
# Stap 5: Call Variants - gatk haplotype caller
#----------------------------------------------

echo "Stap 6: Call Variants - gatk haplotype caller"

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variantcs.vcf

## extract SNPs & INDELS
java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R ${ref} -V ${results}/raw_variantcs.vcf --select-type SNP -O ${results}/raw_snps.vcf

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants -R ${ref} -V ${results}/raw_variantcs.vcf --select-type INDEL -O ${results}/raw_indels.vcf


#--------------------
# Filter Variants - GATK4

## Filter SNPS

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar VariantFiltration \
-R ${ref} \
-V ${results}/raw_snps.vcf \
-O ${results}/filtered_snps.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 60.0" \
-filter-name "MQ_filter" -filter "MQ < 40.0" \
-filter-name "SOR_filter" -filter "SOR > 4.0" \
-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
-genotype-filter-expression "DP < 10" \
-genotype-filter-name "DP_filter" \
-genotype-filter-expression "GQ < 10" \
-genotype-filter-name "GQ_filter"

## Filter INDELS

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar VariantFiltration \
-R ${ref} \
-V ${results}/raw_indels.vcf \
-O ${results}/filtered_indels.vcf \
-filter-name "QD_filter" -filter "QD < 2.0" \
-filter-name "FS_filter" -filter "FS > 200.0" \
-filter-name "MQ_filter" -filter "MQ < 40.0" \
-filter-name "SOR_filter" -filter "SOR > 10.0" \
-genotype-filter-expression "DP < 10" \
-genotype-filter-name "DP_filter" \
-genotype-filter-expression "GQ < 10" \
-genotype-filter-name "GQ_filter"


## Select Variants that PASS filters

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants \
--exclude-filtered \
-V ${results}/filtered_snps.vcf \
-O ${results}/analysis-ready-snps.vcf


java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar SelectVariants \
--exclude-filtered \
-V ${results}/filtered_indels.vcf \
-O ${results}/analysis-ready-indels.vcf
# 
# ## To exclude variants that failed genotype filters
# 
cat ${results}/analysis-ready-snps.vcf | grep -v -E "DP_filter|GQ_filter" > ${results}/analysis-ready-snps-filteredGT.vcf

cat ${results}/analysis-ready-indels.vcf | grep -v -E "DP_filter|GQ_filter" > ${results}/analysis-ready-indels-filteredGT.vcf


# #Stap 6: Annotat Variants - GATK4 Funcotator
# 
# ## Funcotator downloaden 
wget -P ${tools}/ ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator/funcotator_dataSources.v1.7.20200521g.tar.gz
cd ${tools}/
tar -zxf ${tools}/funcotator_dataSources.v1.7.20200521g.tar.gz
cd funcotator_dataSources.v1.7.20200521g
tar -zxf gnomAD_exome.tar.gz
tar -zxf gnomAD_genome.tar.gz
cd
# ## Annotate using Funcotator


java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Funcotator \
--variant ${results}/analysis-ready-snps-filteredGT.vcf \
--reference ${ref} \
--ref-version hg38 \
--data-sources-path ${tools}/funcotator_dataSources.v1.7.20200521g \
--output ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
--output-file-format VCF


java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar Funcotator \
--variant ${results}/analysis-ready-indels-filteredGT.vcf \
--reference ${ref} \
--ref-version hg38 \
--data-sources-path ${tools}/funcotator_dataSources.v1.7.20200521g \
--output ${results}/analysis-ready-indels-filteredGT-functotated.vcf \
--output-file-format VCF


# Extract fields from a VCF file to a tab-delimited table

java -jar ${tools}/gatk-4.5.0.0/gatk-package-4.5.0.0-local.jar VariantsToTable \
-V ${results}/analysis-ready-snps-filteredGT-functotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
-O ${results}/output_snps.table

##
cat ${results}/analysis-ready-snps-filteredGT-functotated.vcf | grep "Funcotation fields are: " | sed 's/|/\t/g' > ${results}/output_curated_variants.txt

cat ${results}/output_snps.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' >> ${results}/output_curated_variants.txt

# Annotate using SnpEff

wget -P ${tools}/ https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
cd ${tools}
unzip snpEff_latest_core.zip
cd snpEff
java -jar snpEff.jar download hg38
cd
java -jar ${tools}/snpEff/snpEff.jar hg38 ${results}/analysis-ready-snps-filteredGT.vcf > ${results}/analysis-ready-snps-filteredGT-annotated.vcf

java -jar ${tools}/snpEff/snpEff.jar hg38 ${results}/analysis-ready-indels-filteredGT.vcf > ${results}/analysis-ready-indels-filteredGT-annotated.vcf

