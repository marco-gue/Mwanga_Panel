# =============================================================================
# 0. Setup
# =============================================================================
library(Biostrings)    # for reverseComplement()

genome_fa <- "genome_blast/NSP306_trifida_chr_v3.fa"    # path to your reference
db_name   <- "~/repos/collaborations/MDP/genome_blast/NSP306_trifida_chr_v3"    # base name for your BLAST DB

# =============================================================================
# 1. (one-time) build BLAST DB if it doesn't exist
# =============================================================================
if (!file.exists(paste0(db_name, ".nin"))) {
  message("Creating BLAST DB: ", db_name)
  system2("makeblastdb",
          args = c("-in",   genome_fa,
                   "-dbtype","nucl",
                   "-out",   db_name))
} else {
  message("BLAST DB already exists: ", db_name)
}

# =============================================================================
# 2. BLAST wrapper (as you provided)
# =============================================================================
blast_my_sweetpotato <- function(query_seq) {
  query_fa <- tempfile(pattern = "qry_", fileext = ".fa")
  writeLines(c(">my_query", query_seq), query_fa)
  
  out_tbl <- tempfile(pattern = "blast_", fileext = ".tsv")
  outfmt_string <- "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
  
  # first pass
  args1 <- c(
    "-task", "blastn-short",
    "-dust", "no",
    "-word_size", "7",
    "-query", query_fa,
    "-db", db_name,
    "-out", out_tbl,
    "-outfmt", shQuote(outfmt_string),
    "-max_target_seqs", "10",
    "-evalue", "1e-5"
  )
  system2("blastn", args = args1, stdout = FALSE, stderr = FALSE)
  
  res <- read.table(out_tbl, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                    col.names = c("query_id","subject_id","perc_identity",
                                  "align_length","q_start","q_end",
                                  "s_start","s_end","evalue","bit_score"))
  return(res)
}

# =============================================================================
# 3. Read primers.txt into a data.frame
# =============================================================================
lines      <- readLines("SSR_A_x_B_pops/primers.txt")
name_idx   <- grep("^PCR_Primers:", lines)
left_idx   <- grep("^\\s*Left:",      lines)
right_idx  <- grep("^\\s*Right:",     lines)

names_vec  <- sub("^PCR_Primers:\\s*", "", lines[name_idx])
left_vec   <- sub("^\\s*Left:\\s*",   "", lines[left_idx])
right_vec  <- sub("^\\s*Right:\\s*",  "", lines[right_idx])

primers_df <- data.frame(
  marker = names_vec,
  fwd    = left_vec,
  rev    = right_vec,
  stringsAsFactors = FALSE
)

# =============================================================================
# 4. Loop over each primer pair, blast fwd+rev-RC, collect results
# =============================================================================
all_hits <- list()

for (i in seq_len(nrow(primers_df))) {
  pr   <- primers_df[i,]
  name <- pr$marker
  
  # forward primer hits
  hits_fwd <- blast_my_sweetpotato(pr$fwd)
  
  # reverse primer RC hits
  rc_seq   <- as.character(reverseComplement(DNAString(pr$rev)))
  hits_rev <- blast_my_sweetpotato(rc_seq)
  
  all_hits[[name]] <- list(
    forward_hits = hits_fwd,
    reverse_hits = hits_rev
  )
}

# =============================================================================
# 5. (Optional) Pairing hits into in-silico amplicons
# =============================================================================
amplicons <- data.frame(
  marker    = character(),
  contig    = character(),
  start     = integer(),
  end       = integer(),
  fwd_evalue= numeric(),
  rev_evalue= numeric(),
  stringsAsFactors = FALSE
)

for (nm in names(all_hits)) {
  fwd <- all_hits[[nm]]$forward_hits
  rev <- all_hits[[nm]]$reverse_hits
  if (nrow(fwd)==0 || nrow(rev)==0) next
  
  # for simplicity, take every combination on same contig,
  # where rev$s_start > fwd$s_start and size plausible
  for (j in seq_len(nrow(fwd))) {
    for (k in seq_len(nrow(rev))) {
      if (fwd$subject_id[j] == rev$subject_id[k]) {
        s1 <- fwd$s_start[j]
        s2 <- rev$s_end[k]
        len <- abs(s2 - s1) + 1
        if (len >= 100 && len <= 400) {
          amplicons <- rbind(amplicons, data.frame(
            marker     = nm,
            contig     = fwd$subject_id[j],
            start      = min(s1,s2),
            end        = max(s1,s2),
            fwd_evalue = fwd$evalue[j],
            rev_evalue = rev$evalue[k],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
}

# View your in-silico SSR amplicons:
print(amplicons)
