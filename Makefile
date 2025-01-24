###############################################################################
# Big Makefile that orchestrates:
#   1) Data download
#   2) Raw alignment
#   3) Variant calling
#   4) Bulk RNA analysis
###############################################################################

# We'll assume dataDownloadScript.sh, raw_data_alignment.sh, variant_calling.sh,
# and bulk_rna.sh are all in the same directory as this Makefile.

.PHONY: all fetch_data alignment variant_calling bulk_rna

# All pipeline steps in sequence
all: fetch_data alignment variant_calling bulk_rna

# 1) Download data
fetch_data:
	chmod +x dataDownloadScript.sh
	./dataDownloadScript.sh

# 2) Run raw data alignment pipeline
alignment:
	chmod +x raw_data_alignment.sh
	./raw_data_alignment.sh

# 3) Run variant calling pipeline
variant_calling:
	chmod +x variant_calling.sh
	./variant_calling.sh

# 4) Run bulk RNA-seq analysis
bulk_rna:
	chmod +x bulk_rna.sh
	./bulk_rna.sh
