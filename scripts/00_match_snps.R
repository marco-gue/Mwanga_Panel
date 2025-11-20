# =============================================================================
# Title:    SNP Mismatch Visualization with Reference and Nucleotide Labels
# Author:   Marcelo Mollinari (adapted)
# Date:     2025-04-22
# Purpose:  Compare multiple sequences against a reference, mark mismatches,
#           display each base, and show the reference atop all others.
# =============================================================================

# ----------------------------
# 1. Load required packages
# ----------------------------
library(tidyverse)

# ----------------------------
# 2. Define reference sequence
# ----------------------------
ref <- 
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTAGTATGATATATATTATAGCCGTGGCAAACAACTCACT"

# ----------------------------
# 3. Split reference into bases
# ----------------------------
ref_bases <- strsplit(ref, "")[[1]]

# ----------------------------
# 4. Create data frame for reference row
# ----------------------------
ref_df <- tibble(
  SeqName  = "Reference",
  Position = seq_along(ref_bases),
  RefBase  = ref_bases,
  Base     = ref_bases,
  Match    = NA
)

# ----------------------------
# 5. Define sequences to compare
# ----------------------------
seqs <- c(
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTACTATGATATATATTATAGCCGTGGCAAACAACTCACT",  # Alt
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTAGTATGATATATATTATAGCCAGCAAACAACTCACTAA",  # RefMatch
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTAGTATGATATATATTATAGCCGCAGCAAACAACTCACT",  # RefMatch
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTAGTATGATATATATTATAGCCGCGGCAAATAACTCACT",  # RefMatch
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTACTATGATATATATTAAAGCCGCGGCAAACAACTCACT",  # AltMatch
  "CCGCTCCCGGACGGTGGCTATTAATTTTGAAATAATAATATGTACTATGATGTATATTATTGCCGCGGCAAACAACTCACT"   # AltMatch
)
names(seqs) <- c("Alt", paste0("Seq_", seq_along(seqs)[-1]))

# ----------------------------
# 6. SNP‐building function
# ----------------------------
get_snp_df <- function(sequence, ref, seqname) {
  ref_split <- strsplit(ref, "")[[1]]
  seq_split <- strsplit(sequence, "")[[1]]
  len <- min(length(ref_split), length(seq_split))
  tibble(
    SeqName  = seqname,
    Position = 1:len,
    RefBase  = ref_split[1:len],
    Base     = seq_split[1:len],
    Match    = (Base == RefBase)
  )
}

# ----------------------------
# 7. Assemble full data frame
# ----------------------------
wanted <- names(seqs)
seq_df <- map2_dfr(seqs[wanted], wanted, ~ get_snp_df(.x, ref, .y))
all_df  <- bind_rows(ref_df, seq_df) %>%
  mutate(SeqName = factor(SeqName, levels = c("Reference", wanted)))

# ----------------------------
# 8. Plot: Reference on top
# ----------------------------
# We'll plot y=SeqName and then set scale_y_discrete so Reference is at the top.
ggplot(all_df, aes(x = Position, y = SeqName, label = Base)) +
  # background tiles
  geom_tile(aes(fill = Match, alpha = is.na(Match)), data = all_df) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 1), guide = "none") +
  scale_fill_manual(
    na.value  = "lightblue",   # reference row
    values    = c("FALSE" = "firebrick", "TRUE" = "gray90"),
    labels    = c("Mismatch", "Match")
  ) +
  # nucleotide text
  geom_text(size = 3, color = "black") +
  # invert the order so Reference is drawn at the top
  scale_y_discrete(limits = rev(levels(all_df$SeqName))) +
  labs(
    title = "Per-base SNP Plot with Reference on Top",
    x     = "Reference Position",
    y     = "Sequence",
    fill  = "Base Match"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid   = element_blank(),
    legend.position = "right"
  )
