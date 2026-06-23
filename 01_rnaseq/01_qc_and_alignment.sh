#!/bin/bash
# RNA-seq QC and STAR alignment pipeline
# Samples: SRR17948859, SRR17948866, SRR17948869, SRR17948871, SRR17948873, SRR17948874
# T2D: SRR17948866, SRR17948874, SRR17948859
# Active controls: SRR17948869, SRR17948871, SRR17948873
set -euo pipefail

RAW_DIR="rnaseq/raw"
QC_DIR="rnaseq/qc"
GENOME_DIR="rnaseq/star_index"
ALIGN_DIR="rnaseq/alignment"
GENOME_FASTA="ref/hg38.fa"
GENOME_GTF="ref/hg38.gtf"
THREADS=8

mkdir -p "$QC_DIR" "$ALIGN_DIR"/{counts,bam,sj,logs}

# --- Step 1: FastQC + MultiQC ---
fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_DIR" -t "$THREADS"
multiqc "$QC_DIR" -o "$QC_DIR" -n multiqc_pre_alignment

# --- Step 2: Build STAR genome index ---
mkdir -p "$GENOME_DIR"
STAR --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$GENOME_FASTA" \
     --sjdbGTFfile "$GENOME_GTF" \
     --sjdbOverhang 100 \
     --runThreadN "$THREADS"

# --- Step 3: STAR alignment ---
for R1 in "$RAW_DIR"/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    STAR --runThreadN "$THREADS" \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --outFileNamePrefix "$ALIGN_DIR/${SAMPLE}_"

    # Organize outputs
    mv "$ALIGN_DIR/${SAMPLE}_Aligned.sortedByCoord.out.bam" "$ALIGN_DIR/bam/"
    mv "$ALIGN_DIR/${SAMPLE}_ReadsPerGene.out.tab" "$ALIGN_DIR/counts/"
    mv "$ALIGN_DIR/${SAMPLE}_SJ.out.tab" "$ALIGN_DIR/sj/"
    mv "$ALIGN_DIR/${SAMPLE}_Log"*.out "$ALIGN_DIR/logs/"

    samtools index "$ALIGN_DIR/bam/${SAMPLE}_Aligned.sortedByCoord.out.bam"
done

echo "Alignment complete. Check logs in $ALIGN_DIR/logs/"
