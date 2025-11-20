# =============================================================================
# Title: Compare Dendrograms of Parents Using G Matrix and Euclidean Distance
# Author: Marcelo Mollinari
# Purpose: Plot two dendrograms side by side — one from genomic relationships
#          and one from Euclidean distances on genotype dosage data.
# =============================================================================

# ----------------------------
# 0. Load Required Libraries
# ----------------------------
library(readr)
library(tidyverse)
library(dendextend)

# ----------------------------
# 1. Load Input Files
# ----------------------------

# Genomic relationship matrix
G.mat <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/G_mat_1328_ind_from_5686_mrks.rds")

# Genotype dosage matrix (markers × samples)
genomat <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_MDP_unique.rds")

# Sample metadata
crosswalk_MDP <- read_csv("~/repos/collaborations/MDP/Rdata_and_spreadsheets/crosswalk_table_sorted_MDP_only.csv")

# ----------------------------
# 2. Get List of Parent IDs
# ----------------------------

moms <- unique(na.omit(crosswalk_MDP$female))
dads <- unique(na.omit(crosswalk_MDP$male))
parents <- unique(c(moms, dads))

# ----------------------------
# 3. Prepare for Plotting
# ----------------------------

# Set layout: 1 row, 2 columns
par(mfrow = c(1, 2), mar = c(7, 4, 4, 2))

# ============================
# A. Dendrogram from G Matrix
# ============================

# Filter parents available in G.mat
parents_G <- parents[parents %in% rownames(G.mat)]
G.parents <- G.mat[parents_G, parents_G]

# Compute distance and clustering
dist_G <- as.dist(1 - G.parents)
hc_G <- hclust(dist_G, method = "average")

# Create dendrogram object
dend_G <- as.dendrogram(hc_G)

# Define labels: Mom or Dad
parent_type_G <- ifelse(parents_G %in% moms, "Mom", "Dad")
names(parent_type_G) <- parents_G
label_colors_G <- ifelse(parent_type_G[labels(dend_G)] == "Mom", "tomato", "steelblue")

# Apply label colors
labels_colors(dend_G) <- label_colors_G
labels_cex(dend_G) <- 1

# Plot
plot(dend_G,
     main = "Dendrogram from G Matrix",
     ylab = "1 - G (Genetic Distance)",
     cex.main = 1.2)
legend("topright", legend = c("Mom", "Dad"), fill = c("tomato", "steelblue"), border = NA)

# ============================
# B. Dendrogram from Euclidean Distance
# ============================

# Filter parents available in genomat
parents_geno <- parents[parents %in% colnames(genomat)]
geno.parents <- t(genomat[, parents_geno])  # rows = individuals

# Compute Euclidean distance and clustering
dist_euc <- dist(geno.parents, method = "euclidean")
hc_euc <- hclust(dist_euc, method = "average")

# Create dendrogram object
dend_euc <- as.dendrogram(hc_euc)

# Define labels: Mom or Dad
parent_type_euc <- ifelse(rownames(geno.parents) %in% moms, "Mom", "Dad")
names(parent_type_euc) <- rownames(geno.parents)
label_colors_euc <- ifelse(parent_type_euc[labels(dend_euc)] == "Mom", "tomato", "steelblue")

# Apply label colors
labels_colors(dend_euc) <- label_colors_euc
labels_cex(dend_euc) <- 1

# Plot
plot(dend_euc,
     main = "Dendrogram from Euclidean Distance",
     ylab = "Euclidean Distance",
     cex.main = 1.2)
legend("topright", legend = c("Mom", "Dad"), fill = c("tomato", "steelblue"), border = NA)
