require(tidyverse)
crosswalk <- read.csv(file = "~/repos/collaborations/MDP/Rdata_and_spreadsheets/crosswalk_table_sorted_MDP_only.csv")
dim(crosswalk)
DT <- data.table::as.data.table(crosswalk)
DT 
# ─── 1. Create family counts from filtered crosswalk ──────────────────────────
family_counts <- crosswalk %>%
  filter(!is.na(female), !is.na(male), female != 0, male != 0) %>%
  group_by(female, male) %>%
  summarise(n.ind = n(), .groups = "drop") %>%
  {
    main <- .
    
    fem_totals <- main %>%
      group_by(female) %>%
      summarise(n.ind = sum(n.ind), .groups = "drop") %>%
      mutate(male = "Total_Half_Sib_2") %>%
      select(female, male, n.ind)
    
    male_totals <- main %>%
      group_by(male) %>%
      summarise(n.ind = sum(n.ind), .groups = "drop") %>%
      mutate(female = "Total_Half_Sib_1") %>%
      select(female, male, n.ind)
    
    grand_total <- tibble(
      female = "Total_Half_Sib_1",
      male = "Total_Half_Sib_2",
      n.ind = sum(main$n.ind, na.rm = TRUE)
    )
    
    bind_rows(main, fem_totals, male_totals, grand_total)
  }

# ─── 2. Define axis orders ────────────────────────────────────────────────────
female_order <- c(
  sort(unique(family_counts$female[family_counts$female != "Total_Half_Sib_1"])),
  "Total_Half_Sib_1"
)

male_order <- c(
  "Total_Half_Sib_2",
  sort(unique(family_counts$male[family_counts$male != "Total_Half_Sib_2"]))
)

# ─── 3. Re-factor and label ───────────────────────────────────────────────────
family_counts <- family_counts %>%
  mutate(
    female = factor(female, levels = female_order),
    male = factor(male, levels = male_order),
    female_lab = fct_relabel(female, ~ str_to_title(str_replace_all(., "_", " "))),
    male_lab = fct_relabel(male, ~ str_to_title(str_replace_all(., "_", " "))),
    is_total = female == "Total_Half_Sib_1" | male == "Total_Half_Sib_2",
    is_grand_total = female == "Total_Half_Sib_1" & male == "Total_Half_Sib_2"
  )

# ─── 4. Label-only data for totals and grand total ────────────────────────────
label_data <- family_counts %>%
  filter(is_total) %>%
  transmute(
    female_lab,
    male_lab,
    label = n.ind
  )

# ─── 5. Main plot data (excluding totals) ─────────────────────────────────────
plot_data <- family_counts %>%
  filter(!is_total)

# ─── 6. Final plot ────────────────────────────────────────────────────────────
ggplot(plot_data, aes(x = female_lab, y = male_lab)) +
  geom_point(aes(size = n.ind, fill = n.ind), shape = 21, colour = "black") +
  geom_text(aes(label = n.ind), colour = "white", size = 3) +
  geom_label(
    data = label_data,
    aes(x = female_lab, y = male_lab, label = label),
    fill = "grey90", fontface = "bold", size = 3.5, label.size = 0.3
  ) +
  scale_fill_gradient(
    name = "n.ind",
    low = "lightblue",
    high = "darkblue",
    na.value = "grey80"
  ) +
  scale_size_area(max_size = 15, guide = "none") +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0),
    axis.text.y = element_text(vjust = 1)
  ) +
  labs(x = NULL, y = NULL)
  