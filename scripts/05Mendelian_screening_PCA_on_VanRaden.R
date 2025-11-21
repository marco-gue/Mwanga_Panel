# =============================================================================
# Title:    PCA on VanRaden Relationship Matrix with Plotly Visualization
# Author:   Marcelo Mollinari
# Purpose:  Compute genomic relationship matrix (G) from dosage matrix, 
#           perform PCA, and visualize results with a 3D interactive plot.
# =============================================================================

# -----------------------
# 0. Load Required Packages
# -----------------------

library(tidyverse)    # Data wrangling and manipulation
library(readr)        # CSV import
library(plotly)       # Interactive 3D plotting
library(viridis)      # Color palettes
library(AGHmatrix)    # Genomic relationship matrices (VanRaden)

# -----------------------
# 1. Import and Preprocess Data
# -----------------------

# Load genotype dosage matrix (rows = markers, columns = samples)
genomat <- readRDS("rdata_and_spreadsheets/dosage.rds")
#u <- apply(genomat, 1, function(x) sum(!is.na(x))/length(x))
#genomat <- genomat[u > 0.8,]
dim(genomat)

# Load sample metadata (e.g. family relationships)
crosswalk_MDP <- read_csv("rdata_and_spreadsheets/crosswalk_table_sorted_only.csv")
crosswalk_MDP <- crosswalk_MDP |> 
  

# Create named vector linking Subject_Barcode to full-sib family
fam.name <- setNames(crosswalk_MDP$full_sib, crosswalk_MDP$Subject_Barcode)
fam.name[is.na(fam.name)] <- names(fam.name[is.na(fam.name)])
fam.name <- fam.name[unique(names(fam.name))]

# Identify parent barcodes
parents <- na.omit(unique(c(crosswalk_MDP$female, crosswalk_MDP$male)))

# -----------------------
# 2. Mendelian Inconsistency Filtering (DISABLED)
# -----------------------

# NOTE: This block is preserved but commented out
#       Use if you need to re-enable Mendelian consistency filtering
ploidy <- 4
geno_classes <- expand_grid(mom = 0:ploidy, dad = 0:ploidy)
valid_genos <- apply(geno_classes, 1, function(z) mappoly::segreg_poly(ploidy, z[1], z[2]) != 0)
geno_labels <- rownames(valid_genos)
seg <- array(FALSE, dim = c(ploidy + 1, ploidy + 1, length(geno_labels)),
             dimnames = list(NULL, NULL, geno_labels))
for (k in seq_along(geno_labels)) {
  seg[,,k] <- matrix(valid_genos[k, ], nrow = ploidy + 1, byrow = TRUE)
}
individuals <- setdiff(colnames(genomat), parents)
for (i in seq_along(individuals)) {
  cat("Processing individual", i, "/", length(individuals), "\n")
  id <- which(crosswalk_MDP$Subject_Barcode == individuals[i])
  mom <- crosswalk_MDP$female[id]
  dad <- crosswalk_MDP$male[id]
  if (any(is.na(c(mom, dad)))) next()
  for (j in 1:nrow(genomat)) {
    mend.dose <- which(seg[genomat[j, mom] + 1, genomat[j, dad] + 1, ]) - 1
    if (!(genomat[j, individuals[i]] %in% mend.dose)) {
      genomat[j, individuals[i]] <- NA
    }
  }
}

genomat[1:5,1:5]

saveRDS(genomat,  "rdata_and_spreadsheets/dosage_barcode_agg_with_Mendelian_filtering.rds")

# Read Mendelian-filtered dosage matrix
genomat <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_MDP_barcode_agg_with_Mendelian_filtering.rds")
genomat <- unique(genomat)
dim(genomat)
saveRDS(genomat,"~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_MDP_unique.rds")


snp.id <- rownames(genomat)
chp <- str_split_fixed(snp.id, "chr|_", 5)[,2:5]
chp[,3] <- sub("^0+", "", chp[,3])
u <- data.frame(SNP = rownames(genomat), 
                Chromosome = as.numeric(chp[,1]), 
                Position = as.numeric(chp[,3]))
CMplot::CMplot(u, type = "p", plot.type = "d", file.output = FALSE)

# -----------------------
# 3. Compute Genomic Relationship Matrix (VanRaden Method)
# -----------------------

G.mat <- Gmatrix(t(genomat), method = "VanRaden", ploidy = 4)

rownames(G.mat) <- setNames(crosswalk_MDP$Subject_Barcode,
                            crosswalk_MDP$Sample_ID)[rownames(G.mat)]
colnames(G.mat) <-setNames(crosswalk_MDP$Subject_Barcode,
                           crosswalk_MDP$Sample_ID)[colnames(G.mat)]

G.mat <- G.mat[!is.na(names(rownames(G.mat))), !is.na(names(colnames(G.mat)))]

saveRDS(G.mat, "rdata_and_spreadsheets/G_mat.rds")

# -----------------------
# 4. Heatmap Plotting Function
# -----------------------

plot_heatmap_by_parent <- function(G, crosswalk, parent_field,
                                   palette_func, palette_args = list(n = 20),
                                   label_angle = 45, label_cex = 0.7) {
  
  ids <- crosswalk$Subject_Barcode[order(crosswalk[[parent_field]])] %>%
    na.omit() %>%
    grep("CIP", ., value = TRUE)
  
  Gsub <- as.matrix(G[ids, ids])
  
  group_vals <- crosswalk[[parent_field]][match(ids, crosswalk$Subject_Barcode)]
  breaks     <- c(which(!duplicated(group_vals)), length(group_vals))
  midpoints  <- (head(breaks, -1) + diff(breaks) / 2) / length(group_vals)
  labels     <- unique(group_vals)
  
  par(mar = c(6, 6, 6, 6), xpd = NA)
  image(Gsub, axes = FALSE, col = do.call(palette_func, palette_args))
  
  usr      <- par("usr")
  x_offset <- usr[1] - 0.02 * diff(usr[1:2])
  y_offset <- usr[3] - 0.02 * diff(usr[3:4])
  
  for (b in breaks) {
    segments(x0 = usr[1], y0 = b / length(group_vals), x1 = usr[2], y1 = b / length(group_vals), col = "darkgray")
    segments(x0 = b / length(group_vals), y0 = usr[3], x1 = b / length(group_vals), y1 = usr[4], col = "darkgray")
  }
  
  axis(1, at = midpoints, labels = FALSE, tck = -0.02)
  text(midpoints, y_offset, labels = labels, srt = label_angle, adj = 1, cex = label_cex)
  
  axis(2, at = midpoints, labels = FALSE, tck = -0.02)
  text(x_offset, midpoints, labels = labels, srt = label_angle, adj = 1, cex = label_cex)
  
  par(xpd = FALSE)
}

# -----------------------
# 5. Generate Heatmaps
# -----------------------

plot_heatmap_by_parent(
  G            = G.mat,
  crosswalk    = crosswalk_MDP,
  parent_field = "female",
  palette_func = viridis::rocket,
  palette_args = list(n = 20, direction = -1)
)

plot_heatmap_by_parent(
  G            = G.mat,
  crosswalk    = crosswalk_MDP,
  parent_field = "male",
  palette_func = viridis::mako,
  palette_args = list(n = 20, direction = -1)
)

# -----------------------
# 5. Generate Heatmaps by Families
# -----------------------

plot_heatmap_by_parent(
  G            = G.mat,
  crosswalk    = crosswalk_MDP,
  parent_field = "female",
  palette_func = viridis::rocket,
  palette_args = list(n = 20, direction = -1)
)

plot_heatmap_by_parent(
  G            = G.mat,
  crosswalk    = crosswalk_MDP,
  parent_field = "male",
  palette_func = viridis::mako,
  palette_args = list(n = 20, direction = -1)
)

# -----------------------
# 6. Perform PCA on G Matrix
# -----------------------

prin_comp <- prcomp(G.mat, scale. = TRUE)
pca_var <- prin_comp$sdev^2
pca_var_percent <- 100 * pca_var / sum(pca_var)

cat("Variance explained by PC1–PC3:", round(sum(pca_var_percent[1:3]), 2), "%\n")

# -----------------------
# 7. Variance Explained Barplot
# -----------------------

barplot(
  height = pca_var_percent[1:15],
  names.arg = paste0("PC", 1:15),
  ylab = "Percentage of Variance Explained",
  xlab = "Principal Components",
  main = "PCA - Variance Explained (Top 15 PCs)",
  col = "steelblue"
)

# -----------------------
# 8. Prepare PCA Results for Plotting
# -----------------------

components <- as.data.frame(prin_comp$x)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

group_labels <- fam.name
gl <- group_labels[rownames(components)]
gl[parents] <- parents

parent.pos <- match(parents, sort(unique(gl)))
pal <- viridis(length(unique(gl)))
pal[parent.pos] <- "red"

# -----------------------
# 9. 3D Interactive PCA Plot (Default Groups)
# -----------------------

fig <- plot_ly(
  data  = components,
  x     = ~PC1, y = ~PC2, z = ~PC3,
  color = ~gl,
  colors= pal,
  text  = ~gl,
  marker= list(size = 5)
) %>%
  layout(
    scene = list(
      bgcolor = "#e5ecf6",
      xaxis   = list(title = "PC1"),
      yaxis   = list(title = "PC2"),
      zaxis   = list(title = "PC3")
    ),
    legend = list(title = list(text = "Group"))
  )

fig

# -----------------------
# 10. Grouping by Parents (Mother / Father)
# -----------------------

components <- as.data.frame(prin_comp$x)
components$PC2 <- -components$PC2
components$PC3 <- -components$PC3

components$Mother <- crosswalk_MDP$female[match(rownames(components), crosswalk_MDP$Subject_Barcode)]
components$Father <- crosswalk_MDP$male[match(rownames(components), crosswalk_MDP$Subject_Barcode)]

is_mom <- rownames(components) %in% unique(crosswalk_MDP$female)
is_dad <- rownames(components) %in% unique(crosswalk_MDP$male)

components$Mother[is_mom] <- paste0("Mom_", rownames(components)[is_mom])
components$Father[is_dad] <- paste0("Dad_", rownames(components)[is_dad])

components$Mother[is.na(components$Mother)] <- "Dad"
components$Father[is.na(components$Father)] <- "Mom"

components$MomParent <- is_mom
components$DadParent <- is_dad

# -----------------------
# 10.1 Plot Grouped by Mother
# -----------------------

mom_levels <- unique(components$Mother)
pal_mom <- c(
  "#000000","#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
  "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2",
  "#f7b6d2","#7f7f7f","#c7c7c7"
)
names(pal_mom) <- c(
  "Dad", "Mom_Huarmeyano", "Huarmeyano", "Mom_Magabali", "Magabali",
  "Mom_Mugande", "Mugande", "Mom_NASPOT_11", "NASPOT_11", "Mom_NASPOT_5",
  "NASPOT_5", "Mom_New_Kawogo", "New_Kawogo", "Mom_Resisto", "Resisto",
  "Mom_Wagabolige", "Wagabolige"
)
pal_mom <- pal_mom[sort(mom_levels)]

fig_mom <- plot_ly(
  data   = components,
  x      = ~PC1, y = ~PC2, z = ~PC3,
  color  = ~Mother,
  colors = pal_mom,
  text   = ~paste("Sample:", rownames(components),
                  "<br>Mom:", Mother,
                  "<br>Dad:", Father),
  marker = list(size = 6),
  visible= "legendonly"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    legend = list(title = list(text = "Mother"))
  )

fig_mom

# -----------------------
# 10.2 Plot Grouped by Father
# -----------------------

dad_levels <- unique(components$Father)
pal_dad <- c(
  "#000000","#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a",
  "#d62728","#ff9896","#9467bd","#c5b0d5","#8c564b","#c49c94","#e377c2",
  "#f7b6d2","#7f7f7f","#c7c7c7"
)
names(pal_dad) <- c(
  "Mom", "Dad_Dimbuka_Bukulula", "Dimbuka_Bukulula", "Dad_Ejumula",
  "Ejumula", "Dad_NASPOT_1", "NASPOT_1", "Dad_NASPOT_10_O", "NASPOT_10_O",
  "Dad_SPK004", "SPK004", "Dad_NASPOT_7", "NASPOT_7", "Dad_NK259L", "NK259L",
  "Dad_NASPOT5/58", "NASPOT5/58"
)

fig_dad <- plot_ly(
  data   = components,
  x      = ~PC1, y = ~PC2, z = ~PC3,
  color  = ~Father,
  colors = pal_dad,
  text   = ~paste("Sample:", rownames(components),
                  "<br>Mom:", Mother,
                  "<br>Dad:", Father),
  marker = list(size = 6),
  visible= "legendonly"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    ),
    legend = list(title = list(text = "Father"))
  )

fig_dad


  
  u <- data.frame(SNP = mrk.id, Chromosome = embedded_to_numeric(x$chrom[mrk.id]), 
                  Position = x$genome.pos[mrk.id])
  CMplot::CMplot(u, type = "p", plot.type = "d", file.output = FALSE)

