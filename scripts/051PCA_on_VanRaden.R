# =============================================================================
# Title:    PCA on VanRaden Relationship Matrix with Heatmap by Parent
# Author:   Marcelo Mollinari
# Purpose:  Compute genomic relationship matrix (G) from dosage data,
#           then plot heatmaps of G sorted by maternal and paternal groups.
# =============================================================================

# 1. Load required libraries -----------------------------------------------
library(tidyverse)    # Data manipulation
library(readr)        # Fast CSV I/O
library(viridis)      # Color palettes (rocket)
library(RColorBrewer) # Spectral palette
library(AGHmatrix)    # VanRaden G matrix construction

# 2. Import and preprocess data --------------------------------------------
# 2.1 Genotype dosage matrix (rows = markers, columns = samples)
dosage_path <- "MDSG_DSp25-10220/Report-DSp25-10220/DSp25-10220_Allele_Dose_Report.csv"
genomat <- read_csv(dosage_path) %>%
  column_to_rownames("MarkerID") %>%
  select(-1:-4)           # drop first 4 metadata columns

# 2.2 Sample metadata (family relationships)
crosswalk_path <- "~/repos/collaborations/MDP/Rdata_and_spreadsheets/crosswalk_table_sorted_MDP_only.csv"
crosswalk <- read_csv(crosswalk_path)

# 2.3 Align column names to Subject_Barcode mapping
barcode_map <- setNames(crosswalk$Subject_Barcode, crosswalk$Sample_ID)
colnames(genomat) <- barcode_map[colnames(genomat)]
# remove any unnamed columns
genomat <- genomat[, !is.na(colnames(genomat))]

# 2.4 Compute VanRaden G matrix (individuals as rows)
ploidy_level <- 6
genomat_t <- t(genomat)
G.mat <- Gmatrix(genomat_t, method = "VanRaden", ploidy = ploidy_level)

# 3. Heatmap plotting function ---------------------------------------------
# Sort by chosen parent field and draw a heatmap with rotated labels
plot_heatmap_by_parent <- function(G, crosswalk, parent_field,
                                   palette_func, palette_args = list(n = 20),
                                   label_angle = 45, label_cex = 0.7) {
  # 3.1 Select IDs and submatrix
  ids <- crosswalk$Subject_Barcode[order(crosswalk[[parent_field]])] %>%
    na.omit() %>%
    grep("UGP", ., value = TRUE)
  Gsub <- as.matrix(G[ids, ids])
  
  # 3.2 Determine group boundaries and labels
  group_vals <- crosswalk[[parent_field]][match(ids, crosswalk$Subject_Barcode)]
  breaks <- c(which(!duplicated(group_vals)), length(group_vals))
  midpoints <- (head(breaks, -1) + diff(breaks)/2) / length(group_vals)
  labels <- unique(group_vals)
  
  # 3.3 Plot heatmap
  par(mar = c(6, 6, 6, 6), xpd = NA)
  image(Gsub, axes = FALSE, col = do.call(palette_func, palette_args))
  
  usr <- par("usr")  # plot region coords
  x_offset <- usr[1] - 0.02 * diff(usr[1:2])
  y_offset <- usr[3] - 0.02 * diff(usr[3:4])
  
  # 3.4 Draw grid lines within plot region
  for (b in breaks) {
    segments(x0 = usr[1], y0 = b/length(group_vals),
             x1 = usr[2], y1 = b/length(group_vals), col = "darkgray")
    segments(x0 = b/length(group_vals), y0 = usr[3],
             x1 = b/length(group_vals), y1 = usr[4], col = "darkgray")
  }
  
  # 3.5 Add rotated axis labels
  # X-axis
  axis(side = 1, at = midpoints, labels = FALSE, tck = -0.02)
  text(x     = midpoints, y = y_offset, labels = labels,
       srt   = label_angle, adj = 1, cex = label_cex)
  # Y-axis
  axis(side = 2, at = midpoints, labels = FALSE, tck = -0.02)
  text(x     = x_offset, y = midpoints, labels = labels,
       srt   = label_angle, adj = 1, cex = label_cex)
  par(xpd = FALSE)
}

# 4. Generate heatmaps -----------------------------------------------------
# 4.1 By mother (female)
plot_heatmap_by_parent(
  G          = G.mat,
  crosswalk  = crosswalk,
  parent_field = "female",
  palette_func = viridis::rocket,
  palette_args = list(n = 20, direction = -1)
)

# 4.2 By father (male)
plot_heatmap_by_parent(
  G            = G.mat,
  crosswalk    = crosswalk,
  parent_field = "male",
  palette_func = viridis::mako,
  palette_args = list(n = 20, direction = -1)
)
