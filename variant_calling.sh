#!/usr/bin/env bash
#
# variant_calling.sh
# 
# 1) bcftools mpileup + call for Hifi subgenomes
# 2) bcftools mpileup + call for RNA-seq
# 3) Filtering, stats, etc.

set -e
set -u

echo "======================================"
echo " Running Variant Calling / VCF Workflow"
echo "======================================"

# 1) Subgenome SNPs with bcftools
bcftools mpileup --threads 4 -Ou --skip-indels \
  -f ./Data/Clover_Data/DNA_Contig1_2.fasta \
  Data_for_VCF_analysis/HIFI_contig_1.bam | \
  bcftools call -mv -Ov -o results/VCF_Files/HIFI_Contig_1.vcf

# ... more steps ...
# E.g. filtering, MultiQC, etc.

echo "Variant calling pipeline completed successfully."

