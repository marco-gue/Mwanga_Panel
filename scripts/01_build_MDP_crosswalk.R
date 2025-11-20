# =============================================================================
# Build and Filter Crosswalk Table for MDP Samples
# Author: Marcelo Mollinari
# =============================================================================

# -----------------------
# 1. Load Required Packages
# -----------------------

library(readxl)
library(dplyr)
library(stringr)

# -----------------------
# 2. Load Input Datasets
# -----------------------

# Load metadata and phenotype/pedigree information
datagrid <- read_excel("~/Mwanga_Panel/plates/DataGrid-24-PEP-02.xlsx")
ped      <- get(load("~/Mwanga_Panel/pedigree/pedigree.rda"))
cip      <- readRDS("~/Mwanga_Panel/phenotype/template_oxa22.rds")
dart     <- read_excel("~/Mwanga_Panel/plates/names_in_dart.xlsx")

# -----------------------
# 3. Construct Base Crosswalk Table
# -----------------------

crosswalk <- datagrid %>%
  transmute(
    Sample_ID       = `Sample ID`,
    Customer_Name   = `Customer Sample Name (optional)`,
    Subject_Barcode = `Subject Barcode (optional)`
  )

# -----------------------
# 4. Merge Pedigree Data by Subject Barcode
# -----------------------

ped_unique <- ped %>%
  distinct(geno, .keep_all = TRUE) %>%
  select(
    Accession_in_Pedigree = geno,
    female,
    male
  )

crosswalk <- crosswalk %>%
  left_join(ped_unique, by = c("Subject_Barcode" = "Accession_in_Pedigree"))

# -----------------------
# 5. Clean and Annotate Crosswalk Table
# -----------------------

crosswalk <- crosswalk %>%
  mutate(
    Synonym_clean = str_extract(Subject_Barcode, "CIP\\d+\\.\\d+"),
    CIP        = Synonym_clean %in% str_extract(cip$geno, "CIP\\d+\\.\\d+"),
    DArT       = Sample_ID %in% dart$MarkerID,
    full_sib      = ifelse(!is.na(female) & !is.na(male),
                           paste0(female, "_x_", male),
                           NA_character_)
  )

# -----------------------
# 6. Finalize and Export Full Crosswalk
# -----------------------

crosswalk <- crosswalk %>%
  select(
    Sample_ID, Customer_Name, Subject_Barcode,
    Synonym_clean,
    female, male, full_sib,
    CIP, DArT
  )

# Preview and save
print(head(crosswalk, 10))
write.csv(
  crosswalk,
  file = "~/Mwanga_Panel/rdata_and_spreadsheets/crosswalk_table_sorted.csv",
  row.names = FALSE
)

# -----------------------
# 7. Filter for Samples Present in DArT and Keep Exclusively MDP Plates 
# -----------------------

# 1. Pull out all non‑zero, non‑NA parents in one go
parents <- crosswalk %>% 
  select(female, male) %>% 
  pivot_longer(everything(), values_drop_na = TRUE, names_to = NULL) %>% 
  distinct(value) %>% 
  pull(value) %>% 
  setdiff(0)

# 2. Filter for in_DArT and keep either parents or full_sib families
crosswalk_filtered <- crosswalk %>%
  filter(DArT,
         Subject_Barcode %in% parents | !is.na(full_sib))
  
# Preview and save filtered table
cat("Filtered crosswalk dimensions:", dim(crosswalk_filtered), "\n")
write.csv(
  crosswalk_filtered,
  file = "~/Mwanga_Panel/rdata_and_spreadsheets/crosswalk_table_sorted_only.csv",
  row.names = FALSE
)
