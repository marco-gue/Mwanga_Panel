# =============================================================================
# sweetpotato_haplo_analysis_parallel.R
#
# Parallelized pipeline to:
#   1) create a BLAST DB from your reference genome (one‑time)
#   2) BLAST each haplotype’s reference sequence
#   3) identify biallelic SNPs and sum read depths for ref/alt alleles
#   4) filter by a depth threshold
#   5) aggregate results into V.ref and V.alt
#
# USAGE:
#   1. Install dependencies:
#        install.packages(c("tibble", "stringr"))
#   2. Make sure NCBI BLAST+ (makeblastdb, blastn) is on your PATH
#   3. Edit the file‐paths below as needed
#   4. Run in R:
#        source("sweetpotato_haplo_analysis_parallel.R")
# =============================================================================

#–– Load libraries
library(parallel)   # for mclapply()
library(stringr)    # for str_split_fixed()
library(tibble)     # for tibble()

# =============================================================================
# 1. Set paths and parameters
# =============================================================================

# Reference genome FASTA
genome_fa <- normalizePath("~/Mwanga_Panel/genome_blast/ATL_v3.asm.fa")

db_name <- "genome_blast/ATL_v3.asm"

# Depth threshold: minimum mean read depth per SNP to keep it
depth.th <- 20

# Number of cores to use (detectCores() - 1 is a good rule of thumb)
ncore <- detectCores() - 5

# =============================================================================
# 2. Create BLAST DB (one‑time step)
# =============================================================================

if (!file.exists(paste0(db_name, ".nin"))) {
  message("Creating BLAST DB: ", db_name)
  system2("makeblastdb",
          args = c("-in", genome_fa,
                   "-dbtype", "nucl",
                   "-out",  db_name))
} else {
  message("BLAST DB already exists: ", db_name)
}


# =============================================================================
# 3. Define BLAST wrapper function
# =============================================================================

blast_my_potato <- function(query_seq) {
  # Write query to FASTA
  query_fa <- tempfile(pattern = "qry_", fileext = ".fa")
  writeLines(c(">my_query", query_seq), query_fa)
  # Format: tabular with selected fields
  outfmt_string <- "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
  # Run blastn
  out_tbl <- tempfile(pattern = "blast_", fileext = ".tsv")
  a1 <-  c(
    "-query", query_fa,
    "-db", db_name,
    "-out", out_tbl,
    "-outfmt", shQuote(outfmt_string),     # Ensure proper quoting of format
    "-max_target_seqs", "10",
    "-evalue", "1e-5"
  ) 
  a2 <-  c(
    "-task", "blastn-short",
    "-dust", "no",
    "-word_size", "7",
    "-query", query_fa,
    "-db", db_name,
    "-out", out_tbl,
    "-outfmt", shQuote(outfmt_string),     # Ensure proper quoting of format
    "-max_target_seqs", "10",
    "-evalue", "1e-5"
  ) 
  
  system2("blastn",
          args = a1,
          stdout = FALSE, 
          stderr = FALSE 
  )
  
  # Read results
  results <- read.table(out_tbl, header = FALSE, sep = "\t",
                        stringsAsFactors = FALSE,
                        col.names = c(
                          "query_id","subject_id",
                          "perc_identity","align_length",
                          "q_start","q_end",
                          "s_start","s_end",
                          "evalue","bit_score"
                        ))
  if(nrow(results) > 0)
    return(results)
  else{
    system2("blastn",
            args = a2,
            stdout = FALSE, 
            stderr = FALSE 
    )
    
    # Read results
    results <- read.table(out_tbl, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE,
                          col.names = c(
                            "query_id","subject_id",
                            "perc_identity","align_length",
                            "q_start","q_end",
                            "s_start","s_end",
                            "evalue","bit_score"
                          ))
    return(results)
  }
}

# Example of use:
# blast_my_sweetpotato("GAGTGTGAAGATTTGGACAAAAGAGGTTAGTTTTTACTGTTATGGCATTTATCTCCTTATAAAATTTTGTATTTTTTTTGT")

# =============================================================================
# 4. Load and prepare input data
# =============================================================================
#–– Read MADC allele count Excel file (only needed once) and save as RDS for faster future loading

# madc_file <- "~/Mwanga_Panel/data/Report-DP24-9865_24-PEP-02/DP24-9865_MADC.csv"
# MADC <- read.csv(madc_file, skip = 7)
# saveRDS(MADC, "~/Mwanga_Panel/rdata_and_spreadsheets/allele_count_file.rds")

# Allele‐count data
MADC <- readRDS(
  "~/Mwanga_Panel/rdata_and_spreadsheets/allele_count_file.rds"
)


# Split by CloneID to get a list of data.frames, one per haplotype
MADC_split_by_microhap <- split(MADC, MADC$CloneID)

# =============================================================================
# 5. Define per‑haplotype processing function
# =============================================================================

process_haplotype <- function(j) {
  x <- MADC_split_by_microhap[[j]]
  # Identify the index of the “reference” allele in this haplotype
  id.ref <- grep("\\|Ref$", x$AlleleID)
  if (length(id.ref) != 1) {
    return(list(ref = NULL, alt = NULL))
  }
  ref.seq <- x$ClusterConsensusSequence[id.ref]
  
  # 1) BLAST the reference sequence
  y <- blast_my_potato(ref.seq)
  
  # If BLAST returned no hits, skip this haplotype
  if (is.null(y) || nrow(y) == 0) {
    return(list(ref = NULL, alt = NULL))
  }
  
  y <- y[which.max(y$perc_identity), , drop = FALSE]
  
  # Safety check in case BLAST returns malformed result
  if (nrow(y) == 0 || is.na(y$s_start) || is.na(y$s_end)) {
    return(list(ref = NULL, alt = NULL))
  }
  
  if(y$s_start > y$s_end){
    temp <- y$s_start
    y$s_start <- y$s_end
    y$s_end <-temp
  }

  # 2) Identify biallelic SNP positions
  sequences <- str_split_fixed(x$AlleleSequence, "", nchar(x$AlleleSequence[1]))
  R_seq    <- str_split(ref.seq, "")[[1]]
  snp.pos  <- which(apply(sequences, 2, function(col) length(unique(col)) == 2))
  if (length(snp.pos) == 0) {
    return(list(ref = NULL, alt = NULL))
  }
  M <- sequences[, snp.pos, drop = FALSE]
  R <- R_seq[snp.pos]
  
  # 3) Extract read-depth matrix (cols 17 onward)
  Z <- x[ , -(1:16), drop = FALSE]
  
  # 4) Split sample indices by allele at each SNP
  W <- apply(M, 2, function(col) split(seq_along(col), col))
  if (!length(W)) {
    return(list(ref = NULL, alt = NULL))
  }
  
  Y.ref <- Y.alt <- NULL
  
  # 5) For each SNP, sum depths for ref & alt groups
  for (i in seq_along(W)) {
    idr <- which(names(W[[i]]) == R[i])
    ida <- which(names(W[[i]]) != R[i])
    if (length(idr) == 1 && length(ida) == 1) {
      ref.allele <- names(W[[i]])[idr]
      alt.allele <- names(W[[i]])[ida]
      pos        <- sprintf("%09d", y$s_start + snp.pos[i] - 1)
      ch         <- y$subject_id
      
      sums_ref <- apply(Z[W[[i]][[idr]], ], 2, sum)
      sums_alt <- apply(Z[W[[i]][[ida]], ], 2, sum)
      
      tpl <- function(sums) {
        tibble(
          snp.id = paste0(ch, "_", pos, "_[", ref.allele, ":", alt.allele, "]_hap_", j),
          !!!as.list(sums)
        )
      }
      
      Y.ref <- rbind(Y.ref, tpl(sums_ref))
      Y.alt <- rbind(Y.alt, tpl(sums_alt))
    }
  }
  
  # 6) Filter by depth threshold
  if (nrow(Y.ref) > 0 && nrow(Y.alt) > 0) {
    keep <- apply(Y.ref[,-1], 1, mean) > depth.th &
      apply(Y.alt[,-1], 1, mean) > depth.th
    Y.ref <- Y.ref[keep, , drop = FALSE]
    Y.alt <- Y.alt[keep, , drop = FALSE]
  }
  
  list(ref = Y.ref, alt = Y.alt)
}

# =============================================================================
# 6. Run in parallel and aggregate results
# =============================================================================

message("Processing ", length(MADC_split_by_microhap), " haplotypes on ", ncore, " cores...")

system.time({
  res_list <- mclapply(
    X        = seq_along(MADC_split_by_microhap),
    FUN      = process_haplotype,
    mc.cores = 1
  )
})

# Combine all per-haplotype data.frames into two global tables
V.ref <- do.call(rbind, lapply(res_list, `[[`, "ref"))
V.alt <- do.call(rbind, lapply(res_list, `[[`, "alt"))

message("Done! 'V.ref' and 'V.alt' are ready with filtered SNP-depth summaries.")

# =============================================================================
# Save results
# =============================================================================
saveRDS(V.ref, "rdata_and_spreadsheets/V_ref.rds")
saveRDS(V.alt, "rdata_and_spreadsheets/V_alt.rds")
# =============================================================================
