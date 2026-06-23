#!/usr/bin/env Rscript
# Integration of RNA-seq DEGs with RRBS differential methylation data.
# Identifies genes that are both differentially expressed AND differentially methylated.

library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(methylKit)
library(poolr)

# ---- Load RNA-seq results ----
rna_seq <- read_tsv("rnaseq/Expression_gene.tsv") %>%
  dplyr::select(gene.id, baseMean, log2FoldChange, padj, hgnc_symbol) %>%
  filter(padj < 0.05) %>%
  rename(Chr = gene.id)

# ---- Load DMP data ----
dmps <- read_tsv("methylation/dmps_all.tsv") %>%
  arrange(chr, start) %>%
  rename(Chr = chr)

# ---- Integration function ----
# For each DEG, find DMPs that fall within the gene body,
# then aggregate methylation changes (Stouffer's method for p-values).
rna_meth_integrating <- function(rna, meth) {
  rna %>%
    mutate(start = map2(gene_start, gene_end,
                        ~meth$start[between(meth$start, .x, .y)])) %>%
    unnest(start) %>%
    left_join(meth, by = c("Chr", "start")) %>%
    drop_na(meth.diff, pvalue) %>%
    group_by(Chr, gene_start, gene_end) %>%
    summarise(
      mean_meth    = mean(meth.diff, na.rm = TRUE),
      liptak_padj  = stouffer(pvalue)$p,
      .groups      = "drop"
    ) %>%
    left_join(rna, ., by = c("Chr", "gene_start", "gene_end"))
}

omics <- rna_meth_integrating(rna_seq, dmps) %>%
  drop_na(mean_meth)

write_tsv(omics, "integration/deg_with_methylation.tsv")

# ---- Summary ----
cat("Genes with both DE and DM:", nrow(omics), "\n")
cat("  Upregulated + hypermethylated:",
    sum(omics$log2FoldChange > 0 & omics$mean_meth > 0), "\n")
cat("  Upregulated + hypomethylated:",
    sum(omics$log2FoldChange > 0 & omics$mean_meth < 0), "\n")
cat("  Downregulated + hypermethylated:",
    sum(omics$log2FoldChange < 0 & omics$mean_meth > 0), "\n")
cat("  Downregulated + hypomethylated:",
    sum(omics$log2FoldChange < 0 & omics$mean_meth < 0), "\n")

cat("Integration complete.\n")
