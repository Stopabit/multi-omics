# Multi-Omics Analysis of Type 2 Diabetes in Skeletal Muscle

**Integrated RNA-seq, RRBS methylation, and ChIP-seq analysis of Type 2 Diabetes (T2D) vs. healthy controls in human vastus lateralis muscle.**

## Biological Background

**Type 2 Diabetes (T2D)** is a metabolic disorder characterized by insulin resistance and impaired glucose homeostasis. Skeletal muscle accounts for ~80% of insulin-stimulated glucose uptake, making it a central tissue in T2D pathophysiology. Molecular changes in muscle — altered gene expression, aberrant DNA methylation, and disrupted chromatin regulation — contribute to the insulin-resistant phenotype, but the interplay between these layers is not fully characterized.

This project integrates three complementary omics technologies to profile the epigenomic and transcriptomic landscape of T2D in **M. vastus lateralis** (quadriceps) biopsies:

| Omics layer | What it measures | Technology | Key question |
|---|---|---|---|
| **RNA-seq** | mRNA abundance | Illumina paired-end sequencing + STAR + DESeq2 | Which genes are differentially expressed in T2D vs. active controls? |
| **RRBS** | CpG methylation at single-base resolution | Reduced Representation Bisulfite Sequencing + Bismark + methylKit | Which CpG sites change methylation in T2D? |
| **ChIP-seq** | Protein–DNA binding (RING1, Polycomb) | Chromatin immunoprecipitation + Bowtie2 + MACS2 | Where does RING1 bind differentially between tissues? |

### Why multi-omics integration?

DNA methylation at promoters typically silences gene expression. If a gene is both **upregulated** (RNA-seq) and **hypomethylated** (RRBS) in T2D, this concordance strengthens the evidence that the methylation change is functionally relevant — not just a passenger event. Conversely, downregulated + hypermethylated genes point to epigenetic silencing. ChIP-seq adds a third layer by identifying Polycomb-mediated chromatin repression, which cooperates with DNA methylation in gene silencing.

## Study Design

**Cohort:** Human subjects from GEO/SRA, M. vastus lateralis biopsies.

| Module | T2D samples | Control samples | Data type |
|---|---|---|---|
| RNA-seq | SRR17948866, SRR17948874, SRR17948859 | SRR17948869, SRR17948871, SRR17948873 | Paired-end mRNA-seq |
| RRBS | SRR17948848, SRR17948856 | SRR17948852, SRR17948854 | Paired-end RRBS |
| ChIP-seq | SRR1055687/88 (BF), SRR1055706/07 (HS) | SRR1055691/92 (BF input), SRR1055710/11 (HS input) | Single-end ChIP-seq |

## Pipeline Overview

```
                    ┌─────────────────────────────────────────┐
                    │           Raw sequencing data            │
                    └────────┬──────────┬──────────┬───────────┘
                             │          │          │
                    ┌────────▼──┐  ┌────▼────┐  ┌──▼────────┐
                    │  RNA-seq   │  │  RRBS   │  │  ChIP-seq │
                    └────────┬──┘  └────┬────┘  └──┬────────┘
                             │          │          │
              FastQC/MultiQC │   Trim Galore      │ FastQC
                             │   (--rrbs)          │
              STAR alignment │   Bismark           │ Bowtie2
              (~95% mapped)  │   (~64% mapped)     │ + dedup
                             │   (~93% BS conv.)   │
              Gene counts    │   Meth. calling     │ CPM BigWig
                             │                     │
              DESeq2         │   methylKit          │ MACS2
              (143 up,       │   (DMPs: Δβ≥10%,    │ peak calling
               90 down)      │    q<0.05)           │ + bdgdiff
                             │                     │
                             │   Genomation         │ HOMER
                             │   annotation         │ annotation
                             │                     │ + motifs
                             └─────────┬───────────┘
                                       │
                              ┌────────▼────────┐
                              │   Integration   │
                              │  DEGs ∩ DMPs    │
                              │  GO enrichment  │
                              └─────────────────┘
```

## Key Results

### 1. RNA-seq — Differential Expression

- **~95% mapping rate** (STAR, hg38)
- **143 upregulated** and **90 downregulated** genes in T2D (|log2FC| > 1.3, padj < 0.05)
- PCA clearly separates T2D from active controls along PC1
- K-means clustering (k=2) recovers the T2D / control groups

### 2. RRBS — DNA Methylation

- **Bisulfite conversion efficiency: ~93%** across all samples
- **Alignment rate: 63–66%** (Bismark, paired-end)
- **Average genome-wide CpG methylation: ~7.7%** (RRBS enriches for CpG-rich regions)
- Bimodal methylation distribution: most CpGs are either fully unmethylated or fully methylated
- **<0.1% CpGs removed** by coverage filtering (lo.count=10, hi.perc=99.9)
- Spearman correlation between all samples: **0.91–0.92** (high concordance)
- PCA separates T2D from controls along PC1
- **DMPs identified** at delta-beta ≥ 10%, Hochberg q < 0.05
- Genomation annotation: **81% of DMPs in non-CpG-island regions**, 11% in shores, 9% in CpG islands
- GO enrichment (GREAT): hyper-DMPs enriched in steroid metabolism and transcription regulation; hypo-DMPs in lipid metabolism

### 3. ChIP-seq — RING1 Binding (Polycomb)

- **RING1** (a core Polycomb Repressive Complex 1 subunit) binding profiled in breast fibroblast (BF) and human skin fibroblast (HS)
- IGV visualization (chr10): clear narrow and broad peaks in experimental samples vs. flat controls
- **MACS2 peak calling:** fragment size d = 118 bp
- **bedtools overlap:** 1 355 common peaks, 2 286 BF-unique, 3 050 HS-unique
- **bdgdiff overlap:** 689 common peaks (more conservative)
- HOMER annotation: peaks distributed across promoters, introns, exons, and intergenic regions
- **Motif discovery** (length 8): tissue-specific motifs identified for BF and HS

### 4. Integration — DEGs with Altered Methylation

RNA-seq DEGs overlapped with RRBS DMPs to identify genes with concordant expression and methylation changes. This highlights candidate genes where epigenetic dysregulation may drive transcriptional changes in T2D muscle.

## Repository Structure

```
.
├── README.md
├── LICENSE
├── .gitignore
│
├── 01_rnaseq/
│   ├── 01_qc_and_alignment.sh         # FastQC → MultiQC → STAR
│   ├── 02_counts_qc.R                 # Coverage distributions, heatmaps, PCA, k-means
│   └── 03_deseq2_deg.R                # DESeq2 differential expression
│
├── 02_rrbs/
│   ├── 01_qc_and_trim.sh              # FastQC → Trim Galore (--rrbs, clip 3bp)
│   └── 02_bismark_alignment.sh        # Bismark alignment + methylation extraction
│
├── 03_differential_methylation/
│   └── 01_methylkit_analysis.R         # methylKit: filtering, clustering, DMPs, genomation, GO
│
├── 04_chipseq/
│   ├── 01_alignment_and_peaks.sh       # Bowtie2 → dedup → BigWig → MACS2 → bedtools
│   ├── 02_annotation_and_motifs.sh     # HOMER annotatePeaks + findMotifs (len 8)
│   └── 03_venn_and_barplots.R          # Venn diagrams + genomic feature barplots
│
└── 05_integration/
    └── 01_rnaseq_methylation_integration.R   # DEG ∩ DMP overlap, Stouffer p-value aggregation
```

## Dependencies

### Command-line tools (conda)

```bash
conda create -n multiomics -c bioconda -c conda-forge \
    fastqc multiqc star samtools bowtie2 \
    trim-galore bismark bedtools macs2 homer
```

### R packages

```r
install.packages(c("ggplot2", "dplyr", "tibble", "readr", "tidyr", "purrr",
                    "factoextra", "gridExtra", "VennDiagram"))

BiocManager::install(c("DESeq2", "biomaRt", "methylKit", "genomation", "rGREAT"))
```

## Reproducing the Analysis

1. Download raw FASTQ files from SRA (accessions listed above)
2. Place reference genome (hg38.fa, hg38.gtf) in `ref/`
3. Run scripts in numerical order within each module:

```bash
# RNA-seq
bash 01_rnaseq/01_qc_and_alignment.sh
Rscript 01_rnaseq/02_counts_qc.R
Rscript 01_rnaseq/03_deseq2_deg.R

# RRBS
bash 02_rrbs/01_qc_and_trim.sh
bash 02_rrbs/02_bismark_alignment.sh

# Differential methylation
Rscript 03_differential_methylation/01_methylkit_analysis.R

# ChIP-seq
bash 04_chipseq/01_alignment_and_peaks.sh
bash 04_chipseq/02_annotation_and_motifs.sh
Rscript 04_chipseq/03_venn_and_barplots.R

# Integration
Rscript 05_integration/01_rnaseq_methylation_integration.R
```

## Methods Summary

| Step | Tool | Key parameters |
|------|------|---------------|
| Read QC | FastQC + MultiQC | — |
| RNA-seq alignment | STAR 2.7+ | `--quantMode GeneCounts`, hg38 |
| Differential expression | DESeq2 | |log2FC| > 1.3, padj < 0.05 |
| RRBS trimming | Trim Galore | `--rrbs --clip_R1 3 --clip_R2 3` |
| Bisulfite alignment | Bismark | Bowtie2 backend, paired-end |
| Methylation calling | Bismark extractor | `--comprehensive --merge_non_CpG` |
| Differential methylation | methylKit | Δβ ≥ 10%, Hochberg q < 0.05 |
| Genomic annotation | genomation | CpG islands + Ensembl gene bodies |
| GO enrichment | rGREAT | basalPlusExt rule |
| ChIP-seq alignment | Bowtie2 | `--very-sensitive`, hg38 |
| Peak calling | MACS2 | `callpeak -f BAM -g hs` |
| Differential peaks | MACS2 bdgdiff + bedtools | intersect / subtract |
| Peak annotation | HOMER annotatePeaks.pl | hg38, `-annStats` |
| Motif discovery | HOMER findMotifs.pl | `-len 8` |
| Multi-omics integration | Custom R | Stouffer's method for DMP p-value aggregation within gene bodies |

## Technologies

STAR | DESeq2 | Bismark | methylKit | MACS2 | HOMER | Bowtie2 | bedtools | genomation | rGREAT | R | Bash

## License

MIT
