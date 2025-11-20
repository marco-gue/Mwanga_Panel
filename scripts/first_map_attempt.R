require(tidyverse)
require(mappoly)

# Genotype dosage matrix: rows = samples, cols = markers
genomat <- readRDS("~/repos/collaborations/MDP/Rdata_and_spreadsheets/dosage_MDP_unique.rds")

# Sample metadata
crosswalk_MDP <- read_csv("~/repos/collaborations/MDP/Rdata_and_spreadsheets/crosswalk_table_sorted_MDP_only.csv")


# your vector
mrk <- rownames(genomat)

# 1. extract the hap‑ID for each marker
#    e.g. from "…_hap_13" grab "13"
hap_id <- sub(".*_hap_(\\d+)$", "\\1", mrk)

# 2. build a named list of markers, one element per haplotype
#    names will be "hap_1", "hap_2", …
hap_list <- split(mrk, paste0("hap_", hap_id))

# 3. group those hap‑lists by how many markers they contain
#    so element "1" has all the hap_… with exactly 1 SNP, "2" has those with 2 SNPs, etc.
nested <- split(hap_list, lengths(hap_list))

# 4. (optional) sort by SNP‐count and give nicer names
counts <- as.integer(names(nested))
nested <- nested[order(counts)]
names(nested) <- paste0(names(nested), "_SNPs")

# inspect
str(nested, max.level=2)
barplot(sapply(nested, length),  
        main = "N. of SNPs per tag", las = 2, col = rev(RColorBrewer::brewer.pal(12, "Spectral")))



# Gather full-sib families
x <- split(crosswalk_MDP$Subject_Barcode, crosswalk_MDP$full_sib)
x <- x[sort(sapply(x, length), decreasing = TRUE)]

i<-1


for(i in 1:length(x)){
  y <- names(x)[i]
  z <- str_split(y, "_x_")
  p1 <- z[[1]][1]
  p2 <- z[[1]][2]
  geno<-genomat[,x[[y]]]
  u <- sapply(strsplit(rownames(geno), "Chr|_"), function(x) x[c(2,3)])
  w<-data.frame(snp_name = rownames(geno),
                P1 = genomat[,p1],
                P2 = genomat[,p2],
                sequence = as.numeric(u[1,]),
                sequence_position = as.numeric(u[2,]),
                geno)
  dat<-mappoly:::table_to_mappoly(w, 6)
  plot(dat)
  dat1 <- filter_missing(dat, type = "individual", filter.thres = .1, inter = FALSE)
  dat1 <- filter_missing(dat1, type = "marker", filter.thres = .1, inter = FALSE)
  dat1 <- filter_individuals(dat1)
  ## 2 SNP per tag
  im <-  intersect(nested$`5_SNPs`[[1]], dat1$mrk.names)
  if(length(im) < 2) next()
  seq.init <- make_seq_mappoly(dat1, im)
  o <- get_genomic_order(seq.init)
  seq.init <- make_seq_mappoly(o)
  tpt <- est_pairwise_rf(input.seq = seq.init)
  m <- rf_list_to_matrix(tpt)
  plot(m, type = "lod")
  bla <- est_rf_hmm_sequential(input.seq = seq.init, 
                               twopt = tpt, start.set = 3,
                               thres.twopt = 20, 
                               thres.hmm = 50)

  ble <- filter_map_at_hmm_thres(bla, 
                                 thres.hmm = 1)
  bli<-reest_rf(ble, tol = 10e-5)
  plot(bli)
  blo <- filter_map_at_hmm_thres(bli, thres.hmm = 0.01)
  plot(blo)
  print(blo, detailed = TRUE)
  blu <- est_full_hmm_with_global_error(blo, error = 0.01, tol = 10e-5)
  plot(blu)
  print(blu, detailed = TRUE)
  z<-calc_genoprob(blu)
  image(t(z$probs[,,3]))
}

u <- est_rf_hmm(seq.init, twopt = tpt, est.given.0.rf = TRUE)
v <- filter_map_at_hmm_thres(u, thres.hmm = 10)
print(v, detailed = TRUE)
w <- est_full_hmm_with_global_error(v, error = 0.01, tol = 10e-5)
print(w, detailed = TRUE)
blu


corrected_matrix <- plot_progeny_dosage_change(list(blu), error=0.01, output_corrected=TRUE) #output corrected
u<-as.numeric(corrected_matrix[,-c(1:4)])
dim(u) <- c(4,22)
u - dat1$geno.dose[rownames(corrected_matrix),]
