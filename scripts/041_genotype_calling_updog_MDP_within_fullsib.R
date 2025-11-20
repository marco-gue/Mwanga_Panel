# =============================================================================
# Title:    Genotype Calling with Updog for MDP Population
# Author:   Marcelo Mollinari
# Date:     2025-04-18
# Purpose:  Process allele-count matrices, perform QC, and call genotypes
#           using the multidog() function from the updog package.
#
# WARNING:  Many full-sib families contain a small number of individuals,
#           which frequently leads to convergence issues during likelihood
#           optimization in Updog. For this reason, this genotype calling
#           approach using the F1 model may not be reliable in its current form.
# =============================================================================

# -----------------------------
# 1. Setup Environment
# -----------------------------
rm(list = ls())                      # Clear workspace
set.seed(12345)                      # Reproducibility
options(stringsAsFactors = FALSE)    # Consistent string handling

library(future)       # Parallelization backend
library(updog)        # Genotype calling
library(readr)        # Fast CSV I/O
library(tidyverse)    # Data wrangling + ggplot2

# Setup parallel processing
ncore <- parallel::detectCores() - 1
plan(multisession, workers = ncore)

# -----------------------------
# 2. Load Input Data
# -----------------------------
# Sample metadata (crosswalk)
crosswalk_MDP <- read_csv(
  "~/repos/collaborations/MDP/Rdata_and_spreadsheets/crosswalk_table_sorted_MDP_only.csv",
  col_types = cols(Sample_ID = col_character())
)

# Allele count matrices
V_ref <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/V_ref.rds")
V_alt <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/V_alt.rds")

# -----------------------------
# 3. Filter to Chromosome SNPs
# -----------------------------
idx_chr <- grepl("Chr", V_ref$snp.id)
stopifnot(all(idx_chr == grepl("Chr", V_alt$snp.id)))

V_ref <- V_ref[idx_chr, ] %>% column_to_rownames("snp.id")
V_alt <- V_alt[idx_chr, ] %>% column_to_rownames("snp.id")

stopifnot(
  all(colnames(V_ref) == colnames(V_alt)),
  all(rownames(V_ref) == rownames(V_alt))
)

# -----------------------------
# 4. Subset to MDP Samples
# -----------------------------
sample_ids <- crosswalk_MDP$Sample_ID
V_ref_MDP <- V_ref[, sample_ids]
V_alt_MDP <- V_alt[, sample_ids]
V_tot_MDP <- V_ref_MDP + V_alt_MDP

# -----------------------------
# 5. Coverage Quality Control
# -----------------------------
mean_counts <- rowMeans(V_tot_MDP)

hist(mean_counts,
     breaks = 200,
     main   = "SNP Coverage Distribution (Mean Total Reads)",
     xlab   = "Mean Total Read Count per SNP")
abline(v = median(mean_counts), col = "red",  lwd = 2)
abline(v = mean(mean_counts),  col = "blue", lwd = 2)

# -----------------------------
# 6. Prepare and Run Updog
# -----------------------------
refmat  <- as.matrix(V_ref_MDP)
sizemat <- as.matrix(V_tot_MDP)

families <- split(crosswalk_MDP$Sample_ID, crosswalk_MDP$full_sib)
family_names <- names(families)
multidog_results <- vector("list", length(families))
names(multidog_results) <- family_names

for (fam in family_names) {
  
  # Extract parent barcodes
  parents <- strsplit(fam, "_x_")[[1]]
  mom_bc <- parents[1]
  dad_bc <- parents[2]
  
  # Get sample IDs
  mom_ids <- crosswalk_MDP$Sample_ID[crosswalk_MDP$Subject_Barcode == mom_bc]
  dad_ids <- crosswalk_MDP$Sample_ID[crosswalk_MDP$Subject_Barcode == dad_bc]
  f1_ids  <- families[[fam]]
  
  # Build ref matrix
  mom_ref <- rowSums(refmat[, mom_ids, drop = FALSE])
  dad_ref <- rowSums(refmat[, dad_ids, drop = FALSE])
  f1_ref  <- refmat[, f1_ids, drop = FALSE]
  colnames(f1_ref) <- crosswalk_MDP$Subject_Barcode[
    match(colnames(f1_ref), crosswalk_MDP$Sample_ID)
  ]
  cur_ref <- cbind(mom_ref, dad_ref, f1_ref)
  colnames(cur_ref)[1:2] <- c(mom_bc, dad_bc)
  
  # Build size matrix
  mom_size <- rowSums(sizemat[, mom_ids, drop = FALSE])
  dad_size <- rowSums(sizemat[, dad_ids, drop = FALSE])
  f1_size  <- sizemat[, f1_ids, drop = FALSE]
  colnames(f1_size) <- crosswalk_MDP$Subject_Barcode[
    match(colnames(f1_size), crosswalk_MDP$Sample_ID)
  ]
  cur_size <- cbind(mom_size, dad_size, f1_size)
  colnames(cur_size)[1:2] <- c(mom_bc, dad_bc)
  
  # Message preview
  message("Running Updog for family: ", fam)
  
  # Try to call genotypes
  mout <- tryCatch({
    multidog(
      refmat  = cur_ref,
      sizemat = cur_size,
      ploidy  = 6,
      p1_id   = mom_bc,
      p2_id   = dad_bc,
      model   = "f1",
      nc      = ncore
    )
  }, error = function(e) {
    warning("multidog() failed for ", fam, ": ", conditionMessage(e))
    return(NULL)
  })
  
  multidog_results[[fam]] <- mout
}

# -----------------------------
# 7. Save Output
# -----------------------------
saveRDS(multidog_results, file = "~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_calling_MDP_f1_model.rds")

# =============================================================================
# End of Script
# =============================================================================
