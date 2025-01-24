#!/usr/bin/env bash
#
# bulk_rna.sh
#
# 1) Merging read counts
# 2) edgeR normalization
# 3) Diff. Expression
# 4) Output generation

set -e
set -u

echo "======================================"
echo " Running Bulk RNA-Seq Pipeline        "
echo "======================================"

# 1) Merging read count files from STAR
samples=$(find results/STAR_output/*_align_contigs_1_2 -name "*ReadsPerGene.out.tab" | sort)
echo "Found samples: $samples"

# R or bash-based merging approach
# e.g. calling an R script: Rscript edgeR_analysis.R ...
# or inline code (fewer lines for example)

# 2) Normalization with edgeR (in R or command-line)
# ...

echo "Bulk RNA analysis completed successfully."
