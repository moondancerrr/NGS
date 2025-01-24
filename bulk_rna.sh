#!/usr/bin/env bash
#
# bulk_rna.sh
#
# This script:
#   1) Finds and merges the STAR "ReadsPerGene.out.tab" files into one matrix
#   2) Invokes R to run edgeR-based bulk RNA analysis (filter, normalize, DGE)
#
# Usage:
#   chmod +x bulk_rna.sh
#   ./bulk_rna.sh
#
# Requirements:
#   - bash, grep, cut, paste
#   - R + packages: edgeR, dplyr, pheatmap, VennDiagram, etc.
#   - The data folder structure from the tutorial, or adapt paths accordingly.

set -e  # Exit on error
set -u  # Treat unset variables as error

echo "============================================================"
echo "        Bulk RNA-seq Analysis: Merge + edgeR DGE           "
echo "============================================================"

###############################################################################
# 1) Merge STAR read counts from "ReadsPerGene.out.tab"
###############################################################################
echo "Step 1) Merging read counts from STAR output..."

# The STAR count files are typically found in:
#   results/STAR_output/<S10_align_contigs_1_2 or TI_align_contigs_1_2>/SampleX.ReadsPerGene.out.tab
# Adjust below if your folder structure is different.
SAMPLES=$(find results/STAR_output/*_align_contigs_1_2 -name "*ReadsPerGene.out.tab" | sort)

if [ -z "$SAMPLES" ]; then
  echo "ERROR: No STAR 'ReadsPerGene.out.tab' files found. Check your paths!"
  exit 1
fi

echo "Found samples:"
echo "$SAMPLES"

# Where to write the merged file
OUT_COUNTS="results/merged_read_counts.tsv"
mkdir -p results

# We'll merge by the second column for each file (unstranded counts).
# The first file sets the gene IDs in col1 and read counts in col2, then we paste subsequent columns.

echo "Creating a temporary merged file..."
FIRST=$(echo "$SAMPLES" | head -1)
# Extract gene IDs and the 2nd col (unstranded) from the first file
paste <(cut -f1 "$FIRST") <(cut -f2 "$FIRST") > "$OUT_COUNTS.tmp"

# For the rest, paste in the 2nd column
for F in $(echo "$SAMPLES" | tail -n +2); do
  paste "$OUT_COUNTS.tmp" <(cut -f2 "$F") > "$OUT_COUNTS.tmp2"
  mv "$OUT_COUNTS.tmp2" "$OUT_COUNTS.tmp"
done

# Build the header line
HEADER="GeneID"
for F in $SAMPLES; do
  BNAME=$(basename "$F" .ReadsPerGene.out.tab)
  HEADER="$HEADER\t$BNAME"
done

# Final merged file
echo -e "$HEADER" > "$OUT_COUNTS"
tail -n +2 "$OUT_COUNTS.tmp" >> "$OUT_COUNTS"
rm -f "$OUT_COUNTS.tmp"

echo "Merged read counts file: $OUT_COUNTS"
echo "------------------------------------------------------------------"

###############################################################################
# 2) Run R for edgeR-based Bulk RNA Analysis
###############################################################################
echo "Step 2) Invoking R to run edgeR analysis..."

# We will read the merged TSV, read a metadata file, then do:
#   - DGEList creation
#   - Filtering low counts
#   - Normalization
#   - Simple design + Contrasts
#   - A test (S10_Treatment vs S10_Control, etc.)
#   - Output some results

Rscript --vanilla << 'EOF'
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(edgeR)
  library(pheatmap)
  library(VennDiagram)
  library(mixOmics)
  library(ggplot2)
})

# ---- Read merged count data ----
counts_file <- "results/merged_read_counts.tsv"
cat("Reading merged counts from:", counts_file, "\n")
Read_counts <- read.delim(counts_file, header=TRUE, row.names=1, check.names=FALSE)

cat("Dimension of merged read_counts:", dim(Read_counts), "\n")
# Example: 366 rows (genes) x 12 columns (samples), or similar.

# ---- Import metadata ----
# The tutorial had a metadata.csv in ../Data/Clover_Data/metadata.csv
# (Condition, Genotype, replicate info).
cat("Reading metadata from ../Data/Clover_Data/metadata.csv\n")
metadata <- read.csv("../Data/Clover_Data/metadata.csv", sep=";", row.names=1, stringsAsFactors=TRUE)
cat("Metadata:\n")
print(metadata)

# Make sure that colnames(Read_counts) match rownames(metadata).
if(!all(colnames(Read_counts) %in% rownames(metadata))){
  cat("WARNING: Mismatch between colnames of Read_counts and rownames of metadata!\n")
  # Ideally reorder or fix. We'll continue but might fail if they differ significantly.
}

# Subset or reorder columns if needed:
# e.g. Read_counts <- Read_counts[, rownames(metadata)]

# ---- Create a DGEList object ----
DGE_Obj <- DGEList(counts=Read_counts, remove.zeros=TRUE)

cat("\nInitial number of genes (including possible zeros):", nrow(DGE_Obj$counts), "\n")

# ---- Filtering lowly expressed genes ----
# The tutorial uses filterByExpr with group=metadata$Group or Condition
# We'll assume metadata has Condition. If not, adapt:
group_vec <- metadata$Condition
keep <- filterByExpr(DGE_Obj, group=group_vec)
DGE_Obj <- DGE_Obj[keep,, keep.lib.sizes=FALSE]
cat("Number of genes after filterByExpr:", nrow(DGE_Obj$counts), "\n")

# ---- Normalization (TMM) ----
DGE_Obj <- calcNormFactors(DGE_Obj, method="TMM")
cat("Norm factors:\n")
print(DGE_Obj$samples$norm.factors)

# ---- Create a design matrix ----
# e.g. Condition is Control/Treatment, 2 levels
# or we might have Genotype+Condition => 4 groups
Condition <- relevel(metadata$Condition, ref="Control") # Set "Control" as baseline
design <- model.matrix(~ Condition, data=metadata)
# If you want Genotype too, do something like:
# design <- model.matrix(~ Genotype + Condition, data=metadata)
cat("\nDesign matrix:\n")
print(design)

# ---- Estimate dispersion, fit QL model ----
DGE_Obj <- estimateDisp(DGE_Obj, design)
fit <- glmQLFit(DGE_Obj, design)

# ---- Quasi-likelihood F-test for Condition Treatment vs Control ----
# The second coefficient in design is ConditionTreatment if the reference is Control
res_qlf <- glmQLFTest(fit, coef=2)
DEG_table <- topTags(res_qlf, n=Inf)$table

# Filter by FDR < 0.05, abs(logFC) > 1
DEG_filtered <- DEG_table[DEG_table$FDR < 0.05 & abs(DEG_table$logFC) > 1,]

cat("\nDifferentially expressed genes:\n")
head(DEG_filtered, n=10)
cat("\nNumber of DE genes passing threshold:", nrow(DEG_filtered), "\n")

write.table(DEG_filtered, file="results/DEG_filtered_condition.txt", quote=FALSE, sep="\t")

# ---- Example: Visualization with a smear plot ----
pdf("results/DEG_smear_plot.pdf")
plotSmear(res_qlf, de.tags=rownames(DEG_filtered))
abline(h=c(-1,1), col="blue")
dev.off()

cat("Done with the basic bulk RNA analysis. Results are in results/DEG_filtered_condition.txt\n")
EOF

echo "------------------------------------------------------------------"
echo "All steps complete. Merged read counts and performed edgeR DGE!"
echo "Results: 'results/DEG_filtered_condition.txt', 'results/DEG_smear_plot.pdf'"
echo "=================================================================="


echo "Bulk RNA analysis completed successfully."
