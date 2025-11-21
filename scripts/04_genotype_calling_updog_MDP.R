# =============================================================================
# Title:    Updog Genotype Calling Pipeline for MDP
# Author:   Marcelo Mollinari
# Date:     2025-04-18
# Purpose:  Load allele‐count data, filter by chromosome, compute total counts,
#           visualize coverage, and run updog’s multidog model.
# =============================================================================

# -----------------------
# 1. Setup environment
# -----------------------
# Clear workspace
rm(list = ls())

# Set a reproducible seed
set.seed(12345)

# Increase string handling precision
options(stringsAsFactors = FALSE)

# Load required packages
library(future)       # for parallel backends
library(updog)        # for genotype calling
library(readr)        # fast CSV import
library(tidyverse)    # dplyr, tibble, ggplot2, etc.

# Detect number of cores to use (leave one free)
ncore <- parallel::detectCores() - 4

# Plan future to use multicore (or multisession on macOS)
plan(multisession, workers = ncore)

# -----------------------
# 2. Import data
# -----------------------
# Crosswalk table mapping MDP samples
crosswalk_MDP <- read_csv(
  "rdata_and_spreadsheets/crosswalk_table_sorted_only.csv",
  col_types = cols(Sample_ID = col_character())
)

# Read in allele‐count matrices for reference and alternate alleles
V_ref <- readRDS("rdata_and_spreadsheets/V_ref.rds")
V_alt <- readRDS("rdata_and_spreadsheets/V_alt.rds")

# -----------------------
# 3. Filter to SNPs on chromosomes
# -----------------------
# Identify rows where snp.id contains "Chr"
idx_ref <- grepl("chr", V_ref$snp.id)
idx_alt <- grepl("chr", V_alt$snp.id)

# Check that filtering indices match
stopifnot(all(idx_ref == idx_alt))

# Subset to only chromosome SNPs
V_ref <- V_ref[idx_ref, ]
V_alt <- V_alt[idx_alt, ]

# Convert SNP IDs to row names for matrix operations
V_ref <- V_ref %>% column_to_rownames("snp.id")
V_alt <- V_alt %>% column_to_rownames("snp.id")

# Verify that dimensions align
stopifnot(
  all(colnames(V_ref) == colnames(V_alt)),
  all(rownames(V_ref) == rownames(V_alt))
)

# -----------------------
# 4. Subset to MDP samples
# -----------------------

fix_names <- function(x) {
  x <- sub("^X", "", x)       # elimina X al inicio
  x <- gsub("\\.", "-", x)    # reemplaza todos los puntos por guiones
  return(x)
}

colnames(V_ref) <- fix_names(colnames(V_ref))
colnames(V_alt) <- fix_names(colnames(V_alt))

# Reorder and subset columns by Sample_ID in crosswalk
samples <- crosswalk_MDP$Sample_ID
V_ref_MDP <- V_ref[, samples]
V_alt_MDP <- V_alt[, samples]

# Compute total read counts per SNP (ref + alt)
V_tot_MDP <- V_ref_MDP + V_alt_MDP

# -----------------------
# 5. Coverage QC plot
# -----------------------
# Mean total count per SNP
q <- apply(V_tot_MDP, 1, mean)
hist(q, breaks = 300)

# Filtering alleles with a mean depth ranging from 10 to 1000:
id <- q > 10 & q < 1000
length(id)
q <- q[id]

# Plot histogram of coverage
hist(q,
     breaks = 300,
     main   = "SNP Coverage Distribution (Mean Total Reads)",
     xlab   = "Mean Total Read Count per SNP")
abline(v = median(q), col = "red", lwd = 2)
abline(v = mean(q),  col = "blue", lwd = 2)

sizemat <- as.matrix(V_tot_MDP[id, ]); dim(sizemat)
refmat <- as.matrix(V_ref_MDP[id, ]); dim(refmat)
altmat <- as.matrix(V_alt_MDP[id, ]); dim(altmat)
save(sizemat, refmat, file = "rdata_and_spreadsheets/refsizemat.rda")
  


# -----------------------
# 6. Prepare matrices for updog
# -----------------------
# updog expects plain matrices
refmat <- as.matrix(V_ref_MDP)
sizemat <- as.matrix(V_tot_MDP)

# -----------------------
# 7. Run multidog model
# -----------------------
mout <- multidog(
  refmat  = refmat,
  sizemat = sizemat,
  ploidy  = 4,
  model   = "norm",
  nc      = ncore
)

# -----------------------
# 8. Visualize results
# -----------------------
# Plot the first 40 SNPs’ genotype fits
plot(mout, indices = 55:60)

saveRDS(mout, file = "rdata_and_spreadsheets/dosage_calling.rds")
genomat <- format_multidog(mout, varname = "geno")
genomat[1:5, 1:5]
saveRDS(genomat, file = "rdata_and_spreadsheets/dosage.rds")

# =============================================================================
# End of script
# =============================================================================
