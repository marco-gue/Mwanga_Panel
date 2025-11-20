# =============================================================================
# Title:    Updog Genotype Calling Pipeline for MDP (Aggregated by Barcode)
# Author:   Marcelo Mollinari
# Date:     2025-04-19
# Purpose:  Load allele‐count data, filter to chromosome SNPs, sum read counts 
#           across all samples sharing the same Subject_Barcode, compute total 
#           counts, visualize coverage, and run updog’s multidog() norm model.
#
# Note:     Summing replicates by barcode increases effective coverage and
#           stabilizes parental/sample estimates.  However, low replication or
#           sparse families may still lead to likelihood convergence issues.
# =============================================================================

# -----------------------
# 1. Setup environment
# -----------------------
rm(list = ls())                      # clear workspace
set.seed(12345)                      # reproducible results
options(stringsAsFactors = FALSE)    # consistent string handling

library(future)       # for parallel backends
library(updog)        # genotype calling
library(readr)        # fast CSV import
library(tidyverse)    # dplyr, tibble, ggplot2, etc.

# detect and reserve one core
ncore <- parallel::detectCores() - 5
plan(multisession, workers = ncore)

# -----------------------
# 2. Import data
# -----------------------
# sample metadata with barcodes
crosswalk_MDP <- read_csv(
  "~/repos/collaborations/MDP/Rdata_and_spreadsheets/crosswalk_table_sorted_MDP_only.csv",
  col_types = cols(
    Sample_ID       = col_character(),
    Subject_Barcode = col_character(),
    full_sib        = col_character()
  )
)

# read allele‐count matrices
V_ref <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/V_ref.rds")
V_alt <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/V_alt.rds")

# -----------------------
# 3. Filter to chromosome SNPs
# -----------------------
chr_idx <- grepl("Chr", V_ref$snp.id)
stopifnot(all(chr_idx == grepl("Chr", V_alt$snp.id)))

V_ref <- V_ref[chr_idx, ] %>% column_to_rownames("snp.id")
V_alt <- V_alt[chr_idx, ] %>% column_to_rownames("snp.id")

stopifnot(
  identical(colnames(V_ref), colnames(V_alt)),
  identical(rownames(V_ref), rownames(V_alt))
)

# -----------------------
# 4. Subset to MDP samples
# -----------------------
sample_ids <- crosswalk_MDP$Sample_ID
V_ref_MDP <- V_ref[, sample_ids]
V_alt_MDP <- V_alt[, sample_ids]

# -----------------------
# 5. Aggregate by barcode
# -----------------------
barcodes <- unique(crosswalk_MDP$Subject_Barcode)

# sum reference counts per SNP across replicates of each barcode
ref_agg <- sapply(barcodes, function(bc) {
  cols <- crosswalk_MDP$Sample_ID[crosswalk_MDP$Subject_Barcode == bc]
  rowSums(V_ref_MDP[, cols, drop = FALSE])
})
rownames(ref_agg) <- rownames(V_ref_MDP)
colnames(ref_agg) <- barcodes

# sum total counts (ref + alt) per SNP across same barcodes
tot_agg <- sapply(barcodes, function(bc) {
  cols <- crosswalk_MDP$Sample_ID[crosswalk_MDP$Subject_Barcode == bc]
  rowSums((V_ref_MDP + V_alt_MDP)[, cols, drop = FALSE])
})
rownames(tot_agg) <- rownames(V_ref_MDP)
colnames(tot_agg) <- barcodes

# Checking read depth increase in the parents
parents <- na.omit(unique(c(crosswalk_MDP$female, crosswalk_MDP$male)))
ref_agg[1:10, parents]
tot_agg[1:10, parents]

# -----------------------
# 6. Coverage QC on aggregated data
# -----------------------
mean_counts_agg <- rowMeans(tot_agg)
hist(mean_counts_agg,
     breaks = 200,
     main   = "Aggregated SNP Coverage (Mean Total Reads)",
     xlab   = "Mean Total Read Count per SNP")
abline(v = median(mean_counts_agg), col = "red",  lwd = 2)
abline(v = mean(mean_counts_agg),  col = "blue", lwd = 2)

# -----------------------
# 7. Prepare matrices for updog
# -----------------------
refmat  <- as.matrix(ref_agg)
sizemat <- as.matrix(tot_agg)

# -----------------------
# 8. Run multidog() norm model
# -----------------------
message("Running multidog() with 'norm' model on aggregated barcodes…")
mout <- tryCatch({
  multidog(
    refmat  = refmat,
    sizemat = sizemat,
    ploidy  = 6,
    model   = "norm",
    nc      = ncore
  )
}, error = function(e) {
  stop("multidog() failed – ", conditionMessage(e))
})

# -----------------------
# 9. Visualize & Save
# -----------------------
# genotype‐fit plots for first 10 SNPs
plot(mout, indices = 1:10)

# save raw multidog output
saveRDS(mout,
        file = "~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_calling_MDP_barcode_agg.rds")
dim(mout$snpdf)
dim(mout$inddf)
head(mout$inddf)

# extract dosage matrix and save
genomat <- format_multidog(mout, varname = "geno")
#idx <- format_multidog(mout, varname = "maxpostprob")
#genomat[idx < 0.5] <- NA
saveRDS(genomat,
        file = "~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_MDP_barcode_agg.rds")
