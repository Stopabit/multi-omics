#!/usr/bin/env Rscript
# Differential methylation analysis with methylKit.
# Input:  Bismark CpG coverage files (*.cov)
# Output: DMPs, per-chromosome plots, genomation annotations, GO enrichment

library(methylKit)
library(genomation)
library(rGREAT)
library(ggplot2)

# ---- Load methylation profiles ----
file.list <- list(
  "rrbs/reports/SRR17948848.bismark.cov",
  "rrbs/reports/SRR17948852.bismark.cov",
  "rrbs/reports/SRR17948854.bismark.cov",
  "rrbs/reports/SRR17948856.bismark.cov"
)

test_data <- methRead(
  file.list,
  sample.id  = list("T2D_1", "Ctrl_1", "Ctrl_2", "T2D_2"),
  assembly   = "hg38",
  treatment  = c(1, 0, 0, 1),
  context    = "CpG",
  mincov     = 10
)

# ---- Coverage and methylation distributions ----
pdf("methylation/meth_coverage_histograms.pdf", width = 6, height = 4)
lapply(test_data, function(.x) {
  getMethylationStats(.x, plot = TRUE, both.strands = FALSE)
  getCoverageStats(.x, plot = TRUE, both.strands = FALSE)
})
dev.off()

# ---- Filter by coverage ----
# lo.count = 10: remove CpGs with < 10 reads
# hi.perc = 99.9: remove top 0.1% (PCR duplication artifacts)
filtered_data <- filterByCoverage(test_data, lo.count = 10, hi.perc = 99.9)

for (i in seq_along(filtered_data)) {
  orig <- nrow(test_data[[i]])
  filt <- nrow(filtered_data[[i]])
  cat(sprintf("%s: %d -> %d CpGs (%.1f%% removed)\n",
              test_data[[i]]@sample.id, orig, filt,
              (orig - filt) / orig * 100))
}

# ---- Merge and assess concordance ----
meth <- unite(filtered_data, destrand = FALSE)

pdf("methylation/sample_correlation.pdf", width = 7, height = 7)
getCorrelation(meth, plot = TRUE, method = "spearman")
dev.off()

pdf("methylation/clustering.pdf", width = 7, height = 7)
clusterSamples(meth, dist = "euclidean")
PCASamples(meth, adj.lim = c(0.4, 0.1), sd.threshold = 0.25,
           filterByQuantile = FALSE)
dev.off()

# ---- Differential methylation ----
diff_meth <- calculateDiffMeth(meth, adjust = "qvalue", mc.cores = 1)

# DMPs: delta-beta >= 10%, Hochberg-corrected q < 0.05
hyper_meth <- getMethylDiff(diff_meth, difference = 10, qvalue = 0.05, type = "hyper")
hypo_meth  <- getMethylDiff(diff_meth, difference = 10, qvalue = 0.05, type = "hypo")
dmps       <- getMethylDiff(diff_meth, difference = 10, qvalue = 0.05)

cat("Hyper-DMPs:", nrow(hyper_meth), "\n")
cat("Hypo-DMPs:", nrow(hypo_meth), "\n")
cat("Total DMPs:", nrow(dmps), "\n")

# ---- Per-chromosome visualization ----
pdf("methylation/diffmeth_per_chr.pdf", width = 10, height = 25)
diffMethPerChr(diff_meth, qvalue.cutoff = 0.05, meth.cutoff = 10, plot = TRUE)
dev.off()

# ---- Genomation annotation ----
cpg_islands <- readBed("ref/cpgIslands.hg38.bed")
gene_annot  <- readBed("ref/hg38.ensembl.fix.bed")

diffCpG_annot <- annotateWithGeneParts(as(dmps, "GRanges"), cpg_islands)
diffGene_annot <- annotateWithGeneParts(as(dmps, "GRanges"), gene_annot)

pdf("methylation/genomation_annotations.pdf", width = 8, height = 5)
plotTargetAnnotation(diffCpG_annot, col = c("red", "blue", "gray"),
                     main = "DMPs in CpG island context")
plotTargetAnnotation(diffGene_annot, col = c("red", "blue", "gray"),
                     main = "DMPs in gene context")
dev.off()

# ---- GO enrichment via GREAT ----
run_great <- function(dmps_subset, label) {
  gr <- as(dmps_subset, "GRanges")
  job <- submitGreatJob(gr, species = "hg38", version = "4",
                        rule = "basalPlusExt",
                        adv_upstream = 3, adv_downstream = 1.5, adv_span = 1)
  tbl <- getEnrichmentTables(job, download_by = "tsv")

  pdf(paste0("methylation/GREAT_", label, ".pdf"), width = 15, height = 5)
  plotRegionGeneAssociationGraphs(job)
  dev.off()

  for (i in seq_along(tbl)) {
    write.table(tbl[[i]],
                file = paste0("methylation/GO_", names(tbl)[i], "_", label, ".tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

run_great(getMethylDiff(hyper_meth), "hyper")
run_great(getMethylDiff(hypo_meth), "hypo")

cat("Differential methylation analysis complete.\n")
