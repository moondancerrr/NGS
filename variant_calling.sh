#!/usr/bin/env bash
#
# variant_calling.sh
#
# This pipeline does the following:
#   1) Subgenome SNP calling from HiFi alignments (Contig1, Contig2, Contig1+2)
#   2) RNA-seq SNP calling for S10 and Tienshan
#   3) Filtering of VCFs
#   4) (Optional) MultiQC or any further analysis
#
# Usage:
#   chmod +x variant_calling.sh
#   ./variant_calling.sh
#
# Requirements:
#   - bcftools >= 1.9
#   - Samtools (for indexing, if needed)
#   - The input .bam files from alignment are assumed to be in:
#       Data_for_VCF_analysis/
#   - The reference genome in ./Data/Clover_Data/DNA_Contig1_2.fasta
#   - Tools on PATH: bcftools, samtools, possibly multiqc

set -e  # Exit script on error
set -u  # Treat unset variables as an error

echo "================================================================"
echo "          Running Variant Calling / VCF Workflow                "
echo "================================================================"

##############################################################################
# Directories and Filenames
##############################################################################
DATA_DIR="./Data_for_VCF_analysis"
REF_DIR="./Data/Clover_Data"
REF_GENOME="${REF_DIR}/DNA_Contig1_2.fasta"
VCF_OUT_DIR="results/VCF_Files"

mkdir -p "$VCF_OUT_DIR"

##############################################################################
# 1) SUBGENOME SNPs: HiFi alignments
#    We call SNPs for:
#      - HIFI_contig_1.bam
#      - HIFI_contig_2.bam
#      - HIFI_contig_1_2.bam
##############################################################################
echo "----------------------------------------------------------------"
echo "1) Subgenome SNP calling using HiFi alignments..."
echo "----------------------------------------------------------------"

# a) HIFI_Contig_1
bcftools mpileup --threads 4 -Ou --skip-indels \
  -f "$REF_GENOME" \
  "$DATA_DIR/HIFI_contig_1.bam" | \
  bcftools call -mv -Ov -o "$VCF_OUT_DIR/HIFI_Contig_1.vcf"

# b) HIFI_Contig_2
bcftools mpileup --threads 4 -Ou --skip-indels \
  -f "$REF_GENOME" \
  "$DATA_DIR/HIFI_contig_2.bam" | \
  bcftools call -mv -Ov -o "$VCF_OUT_DIR/HIFI_Contig_2.vcf"

# c) HIFI_Contig_1_2
bcftools mpileup --threads 4 -Ou --skip-indels \
  -f "$REF_GENOME" \
  "$DATA_DIR/HIFI_contig_1_2.bam" | \
  bcftools call -mv -Ov -o "$VCF_OUT_DIR/HIFI_Contig_1_2.vcf"


##############################################################################
# 2) RNA-SEQ SNPs: S10 and TI alignments
##############################################################################
echo "----------------------------------------------------------------"
echo "2) SNP calling for RNA-seq S10 and Tienshan merges..."
echo "----------------------------------------------------------------"

# a) S10_merged.bam
bcftools mpileup --threads 4 -Ou --skip-indels \
  -f "$REF_GENOME" \
  "$DATA_DIR/RNA_S10_merged.bam" | \
  bcftools call -mv -Ov -o "$VCF_OUT_DIR/RNA_S10_merged.vcf"

# b) TI_merged.bam
bcftools mpileup --threads 4 -Ou --skip-indels \
  -f "$REF_GENOME" \
  "$DATA_DIR/RNA_TI_merged.bam" | \
  bcftools call -mv -Ov -o "$VCF_OUT_DIR/RNA_TI_merged.vcf"


##############################################################################
# 3) Filtering the VCFs
#    Example thresholds: QUAL>150, DP>20, MQ>40
#    Adjust as needed.
##############################################################################
echo "----------------------------------------------------------------"
echo "3) Filtering VCF files with bcftools filter..."
echo "----------------------------------------------------------------"

# Filter each of the five VCFs
for vcf in HIFI_Contig_1.vcf HIFI_Contig_2.vcf HIFI_Contig_1_2.vcf \
           RNA_S10_merged.vcf RNA_TI_merged.vcf
do
  inVCF="$VCF_OUT_DIR/$vcf"
  outVCF="$VCF_OUT_DIR/$(basename $vcf .vcf)_filtered.vcf"

  echo "Filtering $inVCF -> $outVCF ..."
  bcftools filter -i '(QUAL>150) & (DP>20) & (MQ>40)' \
    "$inVCF" -O v > "$outVCF"
done

##############################################################################
# 4) (Optional) Stats, MultiQC, or any further analysis
##############################################################################
# For instance, you could do bcftools stats on each filtered VCF, then gather with multiqc.
# We'll just show an example of bcftools stats. If you want MultiQC to parse them, rename them as .bcfstats etc.

echo "----------------------------------------------------------------"
echo "4) Generating stats for each filtered VCF (example) ..."
echo "----------------------------------------------------------------"

mkdir -p results/vcf_stats

for fvcf in "$VCF_OUT_DIR"/*_filtered.vcf
do
  fname=$(basename "$fvcf" .vcf)
  bcftools stats "$fvcf" > "results/vcf_stats/${fname}.stats"
done

# (Optional) multiqc if you want to aggregate bcftools stats:
# multiqc --outdir results/vcf_stats results/vcf_stats

echo "================================================================"
echo "Variant calling pipeline completed successfully!"
echo "VCFs in: $VCF_OUT_DIR"
echo "Stats in: results/vcf_stats"
echo "================================================================"

