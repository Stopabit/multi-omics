#!/bin/bash
# RRBS QC and trimming pipeline.
# Samples: SRR17948848 (T2D), SRR17948852 (Active), SRR17948854 (Active), SRR17948856 (T2D)
# Tissue: M. Vastus Lateralis
set -euo pipefail

RAW_DIR="rrbs/raw"
QC_PRE="rrbs/qc_pre"
QC_POST="rrbs/qc_post"
TRIM_DIR="rrbs/trimmed"
THREADS=4

mkdir -p "$QC_PRE" "$QC_POST" "$TRIM_DIR"

# --- Step 1: Pre-trimming QC ---
fastqc "$RAW_DIR"/*.fastq.gz -o "$QC_PRE" -t "$THREADS"
multiqc "$QC_PRE" -o "$QC_PRE" -n multiqc_pre_trimming

# --- Step 2: Trim Galore (RRBS mode) ---
# --rrbs: removes filled-in cytosines at fragment ends (MspI digestion artifact)
# --clip_R1/R2 3: removes first 3 non-random nucleotides (visible in per-base content)
# Adapter detection is automatic in Trim Galore
for SAMPLE in SRR17948848 SRR17948852 SRR17948854 SRR17948856; do
    R1=$(ls "$RAW_DIR"/${SAMPLE}*R1*.fastq.gz)
    R2=$(ls "$RAW_DIR"/${SAMPLE}*R2*.fastq.gz)

    trim_galore --rrbs --paired \
                --clip_R1 3 --clip_R2 3 \
                -o "$TRIM_DIR" \
                --basename "$SAMPLE" \
                --gzip "$R1" "$R2"
done

# --- Step 3: Post-trimming QC ---
fastqc "$TRIM_DIR"/*.fq.gz -o "$QC_POST" -t "$THREADS"
multiqc "$QC_POST" -o "$QC_POST" -n multiqc_post_trimming

echo "Trimming complete. Compare $QC_PRE and $QC_POST reports."
