#!/usr/bin/env Rscript
# Differential expression analysis with DESeq2.
# Input:  filtered count matrix + design table
# Output: DEG table, up/downregulated gene lists

library(DESeq2)
library(dplyr)
library(tibble)
library(readr)
library(biomaRt)

# ---- Load data ----
filtered_counts <- read_tsv("rnaseq/filtered_star_counts.tsv") %>%
  column_to_rownames("gene.id") %>% as.data.frame()

design <- tibble(
  Sample = colnames(filtered_counts),
  group  = if_else(Sample %in% c("SRR17948866", "SRR17948874", "SRR17948859"),
                   "T2D", "Active")
)

# ---- Gene annotation via Ensembl ----
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(
  attributes = c("ensembl_gene_id", "ensembl_peptide_id", "hgnc_symbol"),
  filters    = "ensembl_gene_id",
  values     = rownames(filtered_counts),
  mart       = ensembl
)

add_genesym <- function(deseq_res_obj, gene_annot = genemap) {
  idx <- match(rownames(deseq_res_obj), gene_annot$ensembl_gene_id)
  deseq_res_obj$prot        <- gene_annot$ensembl_peptide_id[idx]
  deseq_res_obj$hgnc_symbol <- gene_annot$hgnc_symbol[idx]
  return(deseq_res_obj)
}

# ---- DESeq2 ----
dds <- DESeqDataSetFromMatrix(filtered_counts, design, design = ~ group)
dds <- DESeq(dds)

res <- results(dds) %>%
  as.data.frame() %>%
  add_genesym() %>%
  rownames_to_column("gene.id") %>%
  arrange(log2FoldChange, padj) %>%
  as_tibble()

write_tsv(res, "rnaseq/deseq2_results_full.tsv")

# ---- Filter DEGs ----
# Thresholds: |log2FC| > 1.3, padj < 0.05
dge <- list()
dge$upreg   <- filter(res, log2FoldChange >  1.3 & padj < 0.05)
dge$downreg <- filter(res, log2FoldChange < -1.3 & padj < 0.05)

cat("Upregulated genes (log2FC > 1.3, padj < 0.05):", nrow(dge$upreg), "\n")
cat("Downregulated genes (log2FC < -1.3, padj < 0.05):", nrow(dge$downreg), "\n")

write_tsv(dge$upreg,   "rnaseq/deg_upregulated.tsv")
write_tsv(dge$downreg, "rnaseq/deg_downregulated.tsv")
write_tsv(res, "rnaseq/Expression_gene.tsv")

cat("DEG analysis complete.\n")
