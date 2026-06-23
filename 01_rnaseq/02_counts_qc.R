#!/usr/bin/env Rscript
# Count matrix QC: coverage distributions, sample correlation, PCA, clustering.
# Input:  STAR ReadsPerGene.out.tab files
# Output: filtered count matrix + QC plots

library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(factoextra)

# ---- Load STAR counts ----
ff <- list.files(path = "rnaseq/alignment/counts",
                 pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(ff, read.table, skip = 4)
counts <- as.data.frame(sapply(counts.files, function(x) x[, 2]))
ff <- gsub("[.]ReadsPerGene[.]out[.]tab", "", basename(ff))
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1
counts <- counts[order(row.names(counts)), ]

write_tsv(rownames_to_column(counts, "gene.id"), "rnaseq/table_counts.tsv")

# ---- Experimental design ----
design <- tibble(
  Sample = colnames(counts),
  group  = if_else(Sample %in% c("SRR17948866", "SRR17948874", "SRR17948859"),
                   "T2D", "Active")
)

# ---- Coverage distribution ----
pdf("rnaseq/coverage_distribution.pdf", width = 8, height = 6)
ggplot(data.frame(mean_cov = apply(counts, 1, mean)),
       aes(x = log10(1 + mean_cov))) +
  geom_histogram(color = "black", fill = "lightblue", bins = 50) +
  ggtitle("Distribution of overall mean coverage") +
  labs(x = "log10(mean counts + 1)", y = "Frequency") +
  theme_minimal()

for (i in seq_len(ncol(counts))) {
  p <- ggplot(data.frame(cov = counts[, i]), aes(x = log10(1 + cov))) +
    geom_histogram(color = "black", fill = "lightblue", bins = 50) +
    ggtitle(paste("Coverage distribution:", colnames(counts)[i])) +
    labs(x = "log10(counts + 1)", y = "Frequency") +
    theme_minimal()
  print(p)
}
dev.off()

cat("Genes with zero mean coverage:", sum(apply(counts, 1, mean) == 0), "\n")
cat("Genes with mean coverage >= 10:", sum(apply(counts, 1, mean) >= 10), "\n")

# ---- Filter low-coverage genes ----
filtered_counts <- counts[apply(counts, 1, mean) >= 10, ]
write_tsv(rownames_to_column(filtered_counts, "gene.id"),
          "rnaseq/filtered_star_counts.tsv")

# ---- Sample correlation heatmap ----
cor_logp <- cor(log(1 + counts), method = "pearson")
filtered_cor_logp <- cor(log(1 + filtered_counts), method = "pearson")

pdf("rnaseq/correlation_heatmaps.pdf", height = 6, width = 6)
heatmap(1 - cor_logp, symm = TRUE,
        distfun = function(x) as.dist(x),
        main = "Sample correlation (all genes)")
heatmap(1 - filtered_cor_logp, symm = TRUE,
        distfun = function(x) as.dist(x),
        main = "Sample correlation (filtered)")
dev.off()

# ---- PCA ----
filteredCpm <- sweep(filtered_counts, 2, colSums(filtered_counts), "/") * 1e6
filteredPca <- prcomp(t(filteredCpm), scale = TRUE)
pca_data <- as.data.frame(filteredPca$x[, 1:2]) %>%
  rownames_to_column("Sample") %>%
  left_join(design, by = "Sample")

pct_var <- round(filteredPca$sdev^2 / sum(filteredPca$sdev^2) * 100, 1)

pdf("rnaseq/pca_and_clustering.pdf", height = 6, width = 12)
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample), vjust = -0.5, hjust = 0.5, size = 3) +
  labs(x = paste0("PC1: ", pct_var[1], "%"),
       y = paste0("PC2: ", pct_var[2], "%"),
       colour = "Condition") +
  scale_color_manual(values = c(Active = "cyan3", T2D = "red3")) +
  ggtitle("PCA — T2D vs Active controls") +
  theme_minimal()

# ---- K-means clustering ----
cpm_t <- t(filteredCpm)
km <- kmeans(cpm_t, centers = 2, nstart = 50)
km_plot <- fviz_cluster(km, data = cpm_t) +
  theme_minimal() +
  ggtitle("K-means clustering (k = 2)")

gridExtra::grid.arrange(pca_plot, km_plot, ncol = 2)
dev.off()

cat("QC plots saved to rnaseq/\n")
