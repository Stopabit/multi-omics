#!/bin/bash
# ChIP-seq analysis pipeline: alignment, deduplication, normalization, peak calling.
# Protein: RING1 (Polycomb complex)
# Tissues: BF (breast fibroblast), HS (human skin fibroblast)
# Experiment: SRR1055687/88 (BF ChIP), SRR1055691/92 (BF input)
#             SRR1055706/07 (HS ChIP), SRR1055710/11 (HS input)
set -euo pipefail

RAW_DIR="chipseq/raw"
BAM_DIR="chipseq/bam"
PEAKS_DIR="chipseq/peaks"
BIGWIG_DIR="chipseq/bigwig"
REF="ref/hg38"
CHROM_SIZES="ref/hg38.chrom.sizes"
THREADS=8

mkdir -p "$BAM_DIR" "$PEAKS_DIR" "$BIGWIG_DIR"

# --- Step 1: FastQC ---
fastqc "$RAW_DIR"/*.fastq.gz -t "$THREADS"
multiqc "$RAW_DIR" -o "$RAW_DIR"

# --- Step 2: Bowtie2 alignment + BAM conversion ---
for file in "$RAW_DIR"/*.fq.gz; do
    prefix=$(basename "$file" .fq.gz)
    bowtie2 -p "$THREADS" --no-unal --very-sensitive \
            -x "$REF" -U "$file" \
    | samtools view -@ "$THREADS" -F 4 -Sb \
    | samtools sort -@ "$THREADS" -o "$BAM_DIR/${prefix}.bam"
    samtools index "$BAM_DIR/${prefix}.bam"
done

# --- Step 3: Deduplication ---
for bamfile in "$BAM_DIR"/*.bam; do
    dedup="${bamfile%.bam}.dedup.bam"
    samtools markdup -r "$bamfile" "$dedup"
    samtools index "$dedup"
done

# --- Step 4: CPM normalization + BigWig ---
for bamfile in "$BAM_DIR"/*.dedup.bam; do
    prefix=$(basename "$bamfile" .dedup.bam)
    total_reads=$(samtools view -c "$bamfile")
    scale_factor=$(echo "scale=6; 1000000 / $total_reads" | bc)

    bedtools genomecov -ibam "$bamfile" -bg -scale "$scale_factor" \
        > "$BIGWIG_DIR/${prefix}.bedgraph"
    bedGraphToBigWig "$BIGWIG_DIR/${prefix}.bedgraph" "$CHROM_SIZES" \
        "$BIGWIG_DIR/${prefix}.bw"
    rm "$BIGWIG_DIR/${prefix}.bedgraph"
done

# --- Step 5: MACS2 peak calling ---
macs2 callpeak --bdg \
    --outdir "$PEAKS_DIR" -n RING1_BF -f BAM -g hs \
    -t "$BAM_DIR"/SRR1055687.dedup.bam "$BAM_DIR"/SRR1055688.dedup.bam \
    -c "$BAM_DIR"/SRR1055691.dedup.bam "$BAM_DIR"/SRR1055692.dedup.bam \
    --cutoff-analysis 2> "$PEAKS_DIR/macs2_BF.log"

macs2 callpeak --bdg \
    --outdir "$PEAKS_DIR" -n RING1_HS -f BAM -g hs \
    -t "$BAM_DIR"/SRR1055706.dedup.bam "$BAM_DIR"/SRR1055707.dedup.bam \
    -c "$BAM_DIR"/SRR1055710.dedup.bam "$BAM_DIR"/SRR1055711.dedup.bam \
    --cutoff-analysis 2> "$PEAKS_DIR/macs2_HS.log"

# --- Step 6: Differential peaks (bdgdiff) ---
macs2 bdgdiff \
    --t1 "$PEAKS_DIR/RING1_BF_treat_pileup.bdg" \
    --t2 "$PEAKS_DIR/RING1_HS_treat_pileup.bdg" \
    --c1 "$PEAKS_DIR/RING1_BF_control_lambda.bdg" \
    --c2 "$PEAKS_DIR/RING1_HS_control_lambda.bdg" \
    --outdir "$PEAKS_DIR/bdgdiff" -o-prefix RING1

# --- Step 7: bedtools intersect (shared/unique peaks) ---
bedtools intersect \
    -a "$PEAKS_DIR/RING1_BF_peaks.narrowPeak" \
    -b "$PEAKS_DIR/RING1_HS_peaks.narrowPeak" \
    > "$PEAKS_DIR/common_peaks.bed"

bedtools subtract \
    -a "$PEAKS_DIR/RING1_BF_peaks.narrowPeak" \
    -b "$PEAKS_DIR/RING1_HS_peaks.narrowPeak" \
    > "$PEAKS_DIR/unique_BF_peaks.bed"

bedtools subtract \
    -a "$PEAKS_DIR/RING1_HS_peaks.narrowPeak" \
    -b "$PEAKS_DIR/RING1_BF_peaks.narrowPeak" \
    > "$PEAKS_DIR/unique_HS_peaks.bed"

echo "Peak calling complete."
echo "BF unique: $(wc -l < $PEAKS_DIR/unique_BF_peaks.bed)"
echo "HS unique: $(wc -l < $PEAKS_DIR/unique_HS_peaks.bed)"
echo "Common:    $(wc -l < $PEAKS_DIR/common_peaks.bed)"
