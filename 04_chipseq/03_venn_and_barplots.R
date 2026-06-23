#!/usr/bin/env Rscript
# Venn diagrams and genomic feature barplots for ChIP-seq peak comparison.

library(VennDiagram)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# ---- Venn diagrams: bedtools vs bdgdiff ----
# Peak counts from bedtools intersect
bt_bf_unique <- 2286
bt_hs_unique <- 3050
bt_common    <- 1355

# Peak counts from macs2 bdgdiff
bd_bf_unique <- 931
bd_hs_unique <- 1695
bd_common    <- 689

pdf("chipseq/venn_comparison.pdf", width = 12, height = 6)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

pushViewport(viewport(layout.pos.col = 1))
draw.pairwise.venn(
  area1 = bt_bf_unique + bt_common, area2 = bt_hs_unique + bt_common,
  cross.area = bt_common,
  category = c("RING1_BF", "RING1_HS"),
  fill = c("steelblue", "tomato"), alpha = 0.5,
  cat.cex = 1.2, cex = 1.3
)
grid.text("bedtools intersect", y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport()

pushViewport(viewport(layout.pos.col = 2))
draw.pairwise.venn(
  area1 = bd_bf_unique + bd_common, area2 = bd_hs_unique + bd_common,
  cross.area = bd_common,
  category = c("RING1_BF", "RING1_HS"),
  fill = c("steelblue", "tomato"), alpha = 0.5,
  cat.cex = 1.2, cex = 1.3
)
grid.text("macs2 bdgdiff", y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
popViewport(2)
dev.off()

# ---- Genomic feature barplots from HOMER annStats ----
parse_annstats <- function(path, label) {
  lines <- readLines(path)
  idx <- grep("^Annotation", lines)
  if (length(idx) == 0) return(NULL)
  df <- read.delim(text = lines[idx:length(lines)], header = TRUE, sep = "\t")
  categories <- c("Exon", "Intron", "Promoter", "Intergenic")
  counts <- sapply(categories, function(cat) {
    sum(df$Number.of.peaks[grepl(cat, df$Annotation, ignore.case = TRUE)], na.rm = TRUE)
  })
  tibble(sample = label, feature = categories, count = counts)
}

annstat_files <- list(
  "RING1_BF"    = "chipseq/annotation/RING1_BF_peaks_annStats.txt",
  "RING1_HS"    = "chipseq/annotation/RING1_HS_peaks_annStats.txt",
  "Common"      = "chipseq/annotation/common_peaks_annStats.txt",
  "Unique_BF"   = "chipseq/annotation/unique_BF_peaks_annStats.txt",
  "Unique_HS"   = "chipseq/annotation/unique_HS_peaks_annStats.txt"
)

plot_data <- bind_rows(
  mapply(parse_annstats, annstat_files, names(annstat_files), SIMPLIFY = FALSE)
)

pdf("chipseq/genomic_feature_barplot.pdf", width = 10, height = 6)
ggplot(plot_data, aes(x = sample, y = count, fill = feature)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Number of peaks", fill = "Genomic feature",
       title = "Peak distribution across genomic features") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
dev.off()

cat("Venn diagrams and barplots saved to chipseq/\n")
