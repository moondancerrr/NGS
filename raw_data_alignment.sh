#!/usr/bin/env bash
#
# raw_data_alignment.sh
# 
# 1) FastQC + MultiQC on raw reads
# 2) Reference indexing (samtools faidx, etc.)
# 3) Minimap2 alignment for HiFi reads
# 4) STAR indexing & alignment for RNA-seq
# 5) Qualimap, etc.

set -e  # exit on error
set -u  # treat unset variables as an error

echo "======================================"
echo " Running Raw Data Alignment Workflow "
echo "======================================"

# 1) FASTQC
mkdir -p results/fastqc_output
fastqc -q -o results/fastqc_output ./Data/Clover_Data/*.fastq  > /dev/null 2>&1

# 2) MULTIQC
mkdir -p results/multiqc_output/fastqc_data
multiqc --outdir results/multiqc_output/fastqc_data results/fastqc_output

# 3) Reference indexing
mkdir -p reference_data
cp ./Data/Clover_Data/DNA_Contig1_2.fasta reference_data
samtools faidx reference_data/DNA_Contig1_2.fasta

# 4) Minimapping Hifi reads
mkdir -p results/HIFI_alignment/
minimap2 -a -x map-pb -o results/HIFI_alignment/PacBio_clover_alignment_1_2_mappb.sam \
  reference_data/DNA_Contig1_2.fasta \
  ./Data/Clover_Data/Hifi_reads_white_clover.fastq

##############################################################################
# 5) SORT & INDEX THE SAM FILES -> BAM, THEN QUALIMAP + MULTIQC
##############################################################################
echo "------------------------------------------------------------------"
echo "5) Sorting & indexing the HiFi SAM files..."
echo "------------------------------------------------------------------"

samtools sort results/HIFI_alignment/PacBio_clover_alignment_1_2_mappb.sam \
  -o results/HIFI_alignment/PacBio_clover_alignment_1_2_mappb.sort.bam

samtools sort results/HIFI_alignment/PacBio_clover_alignment_1_2_asm20.sam \
  -o results/HIFI_alignment/PacBio_clover_alignment_1_2_asm20.sort.bam

samtools index results/HIFI_alignment/PacBio_clover_alignment_1_2_mappb.sort.bam
samtools index results/HIFI_alignment/PacBio_clover_alignment_1_2_asm20.sort.bam

echo "------------------------------------------------------------------"
echo "6) Running Qualimap on both sorted BAM files..."
echo "------------------------------------------------------------------"
mkdir -p results/qualimap_output/PacBio_clover_alignment_1_2_mappb
mkdir -p results/qualimap_output/PacBio_clover_alignment_1_2_asm20

qualimap bamqc -bam results/HIFI_alignment/PacBio_clover_alignment_1_2_mappb.sort.bam \
  -outdir results/qualimap_output/PacBio_clover_alignment_1_2_mappb

qualimap bamqc -bam results/HIFI_alignment/PacBio_clover_alignment_1_2_asm20.sort.bam \
  -outdir results/qualimap_output/PacBio_clover_alignment_1_2_asm20

echo "------------------------------------------------------------------"
echo "7) Aggregating Qualimap reports with MultiQC"
echo "------------------------------------------------------------------"
multiqc --outdir results/qualimap_output results/qualimap_output


##############################################################################
# 6) SUBGENOME-SPECIFIC MAPPING (TO CONTIG1 & CONTIG2)
##############################################################################
echo "------------------------------------------------------------------"
echo "8) Mapping HiFi reads separately to Contig1 and Contig2..."
echo "------------------------------------------------------------------"
minimap2 -a -x map-pb -o results/HIFI_alignment/PacBio_clover_alignment_1.sam \
    reference_data/DNA_Contig1.fasta \
    ./Data/Clover_Data/Hifi_reads_white_clover.fastq

minimap2 -a -x map-pb -o results/HIFI_alignment/PacBio_clover_alignment_2.sam \
    reference_data/DNA_Contig2.fasta \
    ./Data/Clover_Data/Hifi_reads_white_clover.fastq

echo "------------------------------------------------------------------"
echo "9) Sorting, indexing, and Qualimap for these subgenome alignments..."
echo "------------------------------------------------------------------"
samtools sort results/HIFI_alignment/PacBio_clover_alignment_1.sam \
  -o results/HIFI_alignment/PacBio_clover_alignment_1.sort.bam
samtools sort results/HIFI_alignment/PacBio_clover_alignment_2.sam \
  -o results/HIFI_alignment/PacBio_clover_alignment_2.sort.bam

samtools index results/HIFI_alignment/PacBio_clover_alignment_1.sort.bam
samtools index results/HIFI_alignment/PacBio_clover_alignment_2.sort.bam

mkdir -p results/qualimap_output/PacBio_clover_alignment_1
mkdir -p results/qualimap_output/PacBio_clover_alignment_2

qualimap bamqc -bam results/HIFI_alignment/PacBio_clover_alignment_1.sort.bam \
  -outdir results/qualimap_output/PacBio_clover_alignment_1

qualimap bamqc -bam results/HIFI_alignment/PacBio_clover_alignment_2.sort.bam \
  -outdir results/qualimap_output/PacBio_clover_alignment_2

echo "------------------------------------------------------------------"
echo "10) MultiQC again for the new Qualimap reports"
echo "------------------------------------------------------------------"
multiqc --outdir results/qualimap_output results/qualimap_output


##############################################################################
# 7) RNA-SEQ MAPPING (STAR): INDEXING, ALIGNMENT, MERGING
##############################################################################
echo "------------------------------------------------------------------"
echo "11) Converting GFF to GTF, then building STAR index on contigs 1&2..."
echo "------------------------------------------------------------------"
gffread -T -o reference_data/white_clover_genes.gtf ../Data/Clover_Data/white_clover_genes.gff

mkdir -p results/STAR_output/indexing_contigs_1_2

# Attempt the first genomeGenerate (may warn about --genomeSAindexNbases)
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir results/STAR_output/indexing_contigs_1_2 \
  --genomeFastaFiles reference_data/DNA_Contig1_2.fasta \
  --sjdbGTFfile reference_data/white_clover_genes.gtf

# Run again with recommended genomeSAindexNbases=9 if warned
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir results/STAR_output/indexing_contigs_1_2 \
  --genomeFastaFiles reference_data/DNA_Contig1_2.fasta \
  --sjdbGTFfile reference_data/white_clover_genes.gtf \
  --genomeSAindexNbases 9


echo "------------------------------------------------------------------"
echo "12) Aligning S10 RNA-seq libraries with STAR..."
echo "------------------------------------------------------------------"
mkdir -p results/STAR_output/S10_align_contigs_1_2

# Loop over S10 R1/R2 fastq pairs
for i in ../Data/Clover_Data/S10*.R1.fastq
do
  PREFIXNAME=$(basename "$i" .R1.fastq)
  echo "-------------------------------------------"
  echo "Aligning PAIRED-END READS: $PREFIXNAME"
  echo "-------------------------------------------"
  STAR --genomeDir results/STAR_output/indexing_contigs_1_2 \
    --runThreadN 8 \
    --runMode alignReads \
    --readFilesIn ../Data/Clover_Data/"$PREFIXNAME".R1.fastq ../Data/Clover_Data/"$PREFIXNAME".R2.fastq \
    --outFileNamePrefix results/STAR_output/S10_align_contigs_1_2/"$PREFIXNAME" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes Standard \
    --quantMode GeneCounts \
    --alignIntronMax 5000
done

echo "------------------------------------------------------------------"
echo "13) Aligning Tienshan RNA-seq libraries with STAR..."
echo "------------------------------------------------------------------"
mkdir -p results/STAR_output/TI_align_contigs_1_2

for i in ../Data/Clover_Data/TI*.R1.fastq
do
  PREFIXNAME=$(basename "$i" .R1.fastq)
  echo "-------------------------------------------"
  echo "Aligning PAIRED-END READS: $PREFIXNAME"
  echo "-------------------------------------------"
  STAR --genomeDir results/STAR_output/indexing_contigs_1_2 \
    --runThreadN 8 \
    --readFilesIn ../Data/Clover_Data/"$PREFIXNAME".R1.fastq ../Data/Clover_Data/"$PREFIXNAME".R2.fastq \
    --outFileNamePrefix results/STAR_output/TI_align_contigs_1_2/"$PREFIXNAME" \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes Standard \
    --quantMode GeneCounts \
    --alignIntronMax 5000
done


##############################################################################
# 8) MULTIQC ON STAR OUTPUT, MERGING BAMs, INDEXING
##############################################################################
echo "------------------------------------------------------------------"
echo "14) MultiQC for S10 and Tienshan STAR alignments"
echo "------------------------------------------------------------------"
mkdir -p results/multiqc_output/S10_STAR_align_1_2
multiqc --outdir results/multiqc_output/S10_STAR_align_1_2 \
  results/STAR_output/S10_align_contigs_1_2

mkdir -p results/multiqc_output/TI_STAR_align_1_2
multiqc --outdir results/multiqc_output/TI_STAR_align_1_2 \
  results/STAR_output/TI_align_contigs_1_2

echo "------------------------------------------------------------------"
echo "15) Merging S10 BAMs -> single sample (samtools merge)"
echo "------------------------------------------------------------------"
mkdir -p results/STAR_output/S10_align_contigs_1_2_merge
samtools merge -f results/STAR_output/S10_align_contigs_1_2_merge/S10.sorted.bam \
  results/STAR_output/S10_align_contigs_1_2/S10_*.sortedByCoord.out.bam

echo "------------------------------------------------------------------"
echo "16) Merging Tienshan BAMs -> single sample"
echo "------------------------------------------------------------------"
mkdir -p results/STAR_output/TI_align_contigs_1_2_merge
samtools merge -f results/STAR_output/TI_align_contigs_1_2_merge/TI.sorted.bam \
  results/STAR_output/TI_align_contigs_1_2/TI_*.sortedByCoord.out.bam

echo "------------------------------------------------------------------"
echo "17) Indexing merged BAMs"
echo "------------------------------------------------------------------"
samtools index results/STAR_output/S10_align_contigs_1_2_merge/S10.sorted.bam
samtools index results/STAR_output/TI_align_contigs_1_2_merge/TI.sorted.bam


echo "Raw data alignment pipeline completed successfully."

