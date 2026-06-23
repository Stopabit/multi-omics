#!/bin/bash
# Peak annotation with HOMER and motif discovery.
set -euo pipefail

PEAKS_DIR="chipseq/peaks"
ANNOT_DIR="chipseq/annotation"
MOTIF_DIR="chipseq/motifs"

mkdir -p "$ANNOT_DIR" "$MOTIF_DIR"

# --- Step 1: Annotate peaks with HOMER ---
for PEAK_FILE in \
    "$PEAKS_DIR/RING1_BF_peaks.narrowPeak" \
    "$PEAKS_DIR/RING1_HS_peaks.narrowPeak" \
    "$PEAKS_DIR/common_peaks.bed" \
    "$PEAKS_DIR/unique_BF_peaks.bed" \
    "$PEAKS_DIR/unique_HS_peaks.bed"; do

    NAME=$(basename "$PEAK_FILE" | sed 's/\.\(narrowPeak\|bed\)$//')
    annotatePeaks.pl "$PEAK_FILE" hg38 \
        -annStats "$ANNOT_DIR/${NAME}_annStats.txt" \
        > "$ANNOT_DIR/${NAME}_annotated.tsv"
done

# --- Step 2: Motif discovery (length 8) ---
# Extract gene lists (column 13 from HOMER annotation)
for TISSUE in BF HS; do
    cut -f16 "$ANNOT_DIR/RING1_${TISSUE}_peaks_annotated.tsv" \
        | tail -n +2 | sort -u \
        > "$MOTIF_DIR/genes_${TISSUE}.txt"

    findMotifs.pl "$MOTIF_DIR/genes_${TISSUE}.txt" human \
        "$MOTIF_DIR/motifs_${TISSUE}" -len 8
done

echo "Annotation and motif analysis complete."
echo "Compare motif results in $MOTIF_DIR/motifs_BF/ and $MOTIF_DIR/motifs_HS/"
