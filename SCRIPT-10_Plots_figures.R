### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load required libraries.
###############################
library(lubridate)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(ggrepel)

### Load variant data.
#########################
Collated_variant <- read.csv("Results_variant_recoder_protein_df_final.csv")

### Load the updated Oxford classifications from script 5. 
#############################################################
Oxford_classifications <- read.csv("Oxford_classifications_for_all_lab_variants - EDIT.csv")

# Apply permanent patch to variant c.388G>A for Oxford in the collated cohort
Oxford_classifications <- Oxford_classifications %>%
  mutate(
    Final_classification = if_else(
      grepl("c\\.388G>A$", Variant_cnomen),
      "LP",
      Final_classification
    )
  )

## Save the updated reference file version. 
write.csv(Oxford_classifications, "Oxford_classifications_for_all_lab_variants - EDIT2.csv", row.names = TRUE, quote = TRUE)

### Load the functional data.
################################
Variants_combined_with_APC_functional_data_unique <- read.csv("Variants_combined_with_APC_functional_data_unique.csv")

# Apply patch to variant c.388G>A for Oxford in the functional dataset
Variants_combined_with_APC_functional_data_unique <- Variants_combined_with_APC_functional_data_unique %>%
  mutate(
    Lab.variant.classification = if_else(
      Lab == "Oxford" & grepl("c\\.388G>A$", Final_cnomen),
      "Likely pathogenic",
      Lab.variant.classification
    )
  )

### Load all individual criteria data frames.
#################################################
Functional_evidence_result <- read.csv("Functional_evidence_result_KCNH2_PS3_BS3.csv")
Allele_freq_criterion_result <- read.csv("Allele_freq_criterion_result_GnomAD_260426.csv")
REVEL_evidence_result <- read.csv("REVEL_evidence_result_KCNH2.csv")
SpliceAI_evidence_result <- read.csv("SpliceAI_evidence_result_KCNH2.csv")
PM1_domain_criterion_results <- read.csv("KCNH2_PM1_hotspot_evaluations.csv")


### CRITICAL UPSTREAM ANONYMIZATION ENGINE
##########################################
# Define your exact laboratory translation dictionary mapping profiles
# This has been edited to keep the labs anon for the repository. This will not run without editing. 
lab_lookup <- c(
  "A"   = "Lab 1",
  "B"    = "Lab 2",
  "C"   = "Lab 3",
  "D" = "Lab 4",
  "E"     = "Lab 5"
)

# FIXED: Native character formatting handles label lookups reliably across packages
Collated_variant$Lab <- as.character(lab_lookup[Collated_variant$Lab])
Variants_combined_with_APC_functional_data_unique$Lab <- as.character(lab_lookup[Variants_combined_with_APC_functional_data_unique$Lab])


### Variant distribution.
#########################
## Create a temporary clean nomenclature column for precise cross-referencing
Collated_variant$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", Collated_variant$Final_cnomen)

## Left join the relevant classification and date columns from the Oxford file
Collated_variant <- Collated_variant %>%
  left_join(Oxford_classifications %>% select(Variant_cnomen, Final_classification, Classified_date), 
            by = c("Match_cnomen" = "Variant_cnomen"))

## Standardise the classification terminology layout style to match other labs
Collated_variant <- Collated_variant %>%
  mutate(
    Final_classification = case_when(
      Final_classification == "P"  ~ "Pathogenic",
      Final_classification == "LP" ~ "Likely pathogenic",
      Final_classification == "LB" ~ "Likely benign",
      Final_classification == "B"  ~ "Benign",
      TRUE                          ~ Final_classification  # Retains VUS, Risk allele, or unchanged metrics
    )
  )

## Identify rows belonging to the Oxford lab (now "Lab 5") where a valid classification exists
Oxford_rows <- Collated_variant$Lab == "Lab 5" & 
  !(Collated_variant$Final_classification == "N/A") & 
  !(Collated_variant$Final_classification == "N/A.") & 
  !is.na(Collated_variant$Final_classification) &
  Collated_variant$Final_classification != ""

## Update the variant classification and classification date columns for the valid Oxford rows
Collated_variant$Lab.variant.classification[Oxford_rows] <- Collated_variant$Final_classification[Oxford_rows]
Collated_variant$Classification.date[Oxford_rows]        <- Collated_variant$Classified_date[Oxford_rows]

## Clean up the dataframe by removing the temporary and joined lookup columns
Collated_variant$Match_cnomen <- NULL
Collated_variant$Final_classification <- NULL
Collated_variant$Classified_date <- NULL

# Apply permanent patch to variant c.388G>A for Oxford in the collated cohort
Collated_variant <- Collated_variant %>%
  mutate(
    Lab.variant.classification = if_else(
      Lab == "Lab 5" & grepl("c\\.388G>A$", Final_cnomen),
      "Likely pathogenic",
      Lab.variant.classification
    )
  )

### Edit collated dataset to account for duplicates in Lab X data. 
##################################################################
# Ensure the temporary matching column exists for grouping
Collated_variant$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", Collated_variant$Final_cnomen)

# 1. Separate non-Brompton data to ensure it remains completely unaffected
# Edited for the repository, won't run unless fixed. 
non_brompton_df <- Collated_variant %>% filter(Lab != "Lab X")

# 2. Isolate Brompton data and parse the date and lab number fields for numeric comparison
brompton_df <- Collated_variant %>% 
  filter(Lab == "Lab 3") %>%
  mutate(
    Parsed_Date = parse_date_time(Classification.date, orders = c("dmy", "ymd", "mdy")),
    Numeric_Lab_Num = as.numeric(gsub("[^0-9]", "", Lab.number))
  ) %>%
  mutate(
    Numeric_Lab_Num = if_else(is.na(Numeric_Lab_Num), 0, Numeric_Lab_Num)
  )

# 3. Clean Brompton data sequentially
brompton_cleaned <- brompton_df %>%
  group_by(Match_cnomen) %>%
  mutate(
    Any_Date_Missing = any(is.na(Classification.date) | trimws(Classification.date) == ""),
    Max_Group_Date = if(all(is.na(Parsed_Date))) as.POSIXct(NA) else max(Parsed_Date, na.rm = TRUE)
  ) %>%
  filter(
    Any_Date_Missing | (!is.na(Parsed_Date) & Parsed_Date == Max_Group_Date)
  ) %>%
  mutate(
    Max_Lab_Num = max(Numeric_Lab_Num, na.rm = TRUE)
  ) %>%
  filter(
    Numeric_Lab_Num == Max_Lab_Num
  ) %>%
  slice(1) %>%
  ungroup() %>%
  select(-Parsed_Date, -Numeric_Lab_Num, -Any_Date_Missing, -Max_Group_Date, -Max_Lab_Num)

# 4. Recombine
Collated_variant <- rbind(non_brompton_df, brompton_cleaned)
Collated_variant$Match_cnomen <- NULL

## Save updated variant list. 
write.csv(Collated_variant, "Results_variant_recoder_protein_df_plots.csv", row.names = TRUE, quote = TRUE)

### Tidy Visualisations Setup.
##############################
df_plots <- read.csv("Results_variant_recoder_protein_df_plots.csv", stringsAsFactors = FALSE)

# Generate uniform clean dataset base configuration
df_base_cleaned <- df_plots %>%
  filter(!is.na(Lab) & Lab != "") %>%
  mutate(
    Classification = trimws(Lab.variant.classification),
    Classification = case_when(
      is.na(Classification) | Classification == "" | Classification == "NaN" ~ "Classification not provided",
      Classification %in% c("VUS", "Variant of uncertain significance", "Uncertain Significance") ~ "VUS",
      Classification %in% c("Likely Pathogenic", "Likely pathogenic", "LP") ~ "Likely Pathogenic",
      Classification %in% c("Pathogenic", "pathogenic", "P") ~ "Pathogenic",
      Classification %in% c("Likely Benign", "Likely benign", "LB") ~ "Likely Benign",
      Classification %in% c("Benign", "benign", "B") ~ "Benign",
      Classification == "Risk allele" ~ "Risk Allele",
      TRUE ~ Classification
    )
  )

# Generate uniform clean dataset missense configuration
df_missense_cleaned <- df_base_cleaned %>%
  filter(
    grepl("^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$", Final_pnomen) &
      !grepl("^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)", Final_pnomen) &
      grepl("^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$", Final_cnomen)
  )

## Number of variants from each lab and proportions of classifications. 
#######################################################################
df_fig1_data <- df_base_cleaned %>%
  distinct(Lab, Final_cnomen, .keep_all = TRUE) 

df_totals <- df_fig1_data %>%
  group_by(Lab) %>%
  summarise(Total_Count = n(), .groups = 'drop')

clinical_order <- c("Benign", "Likely Benign", "VUS", "Risk Allele", "Likely Pathogenic", "Pathogenic", "Classification not provided")
available_classes <- unique(df_fig1_data$Classification)
final_order <- c(clinical_order[clinical_order %in% available_classes], 
                 available_classes[!available_classes %in% clinical_order])

df_fig1_data$Classification <- factor(df_fig1_data$Classification, levels = final_order)

# Generate Figure 1
p1 <- ggplot(df_fig1_data, aes(x = Lab, fill = Classification)) +
  geom_bar(position = "stack", color = "black", width = 0.7) +
  geom_text(data = df_totals, aes(x = Lab, y = Total_Count, label = Total_Count),
            vjust = -0.5, size = 4.5, fontface = "bold", color = "black", inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) + 
  scale_fill_manual(
    values = c(
      "Benign" = "#2c7bb6",
      "Likely Benign" = "#abd9e9",
      "VUS" = "#ffffbf",
      "Risk Allele" = "#b2abd2",
      "Likely Pathogenic" = "#fdae61",
      "Pathogenic" = "#d7191c",
      "Classification not provided" = "#969696"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Unique KCNH2 Variants by Submitting Laboratory",
    x = "Submitting Laboratory",
    y = "Number of Unique Variants",
    fill = "Laboratory Classification"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", color = "black", margin = margin(t = 6)),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(color = "black", size = 0.5), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

print(p1)
ggsave("Figure_01_Laboratory_Variant_Distribution.png", plot = p1, width = 10, height = 6, dpi = 300)

### Remove rows with empty classifications. 
###########################################
# Manchester filter targets "Lab 4" cleanly
df_missense_cleaned <- df_missense_cleaned %>%
  filter(!(Lab == "Lab 4" & (is.na(Lab.variant.classification) | trimws(Lab.variant.classification) == "")))

## Overall variant number grouped by classification. 
####################################################
df_unique_variants <- df_missense_cleaned %>%
  group_by(Final_cnomen) %>%
  summarise(
    valid_classes = list(unique(na.omit(Classification))),
    num_classes = length(unlist(valid_classes)),
    .groups = 'drop'
  ) %>%
  mutate(
    Consensus_Classification = case_when(
      num_classes > 1 ~ "Conflicting Classification",
      num_classes == 1 ~ sapply(valid_classes, function(x) x[1]),
      TRUE ~ "Classification not provided"
    )
  ) %>%
  filter(!Consensus_Classification %in% c("Classification not provided", "Risk Allele"))

axis_order <- c("Benign", "Likely Benign", "VUS", "Conflicting Classification", "Likely Pathogenic", "Pathogenic")
df_unique_variants$Consensus_Classification <- factor(df_unique_variants$Consensus_Classification, levels = axis_order)

p2 <- ggplot(df_unique_variants, aes(x = Consensus_Classification, fill = Consensus_Classification)) +
  geom_bar(color = "black", width = 0.6, show.legend = FALSE) + 
  geom_text(stat = 'count', aes(label = ..count..), 
            vjust = -0.5, size = 4.5, fontface = "bold", color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  scale_fill_manual(
    values = c(
      "Benign" = "#2c7bb6",
      "Likely Benign" = "#abd9e9",
      "VUS" = "#ffffbf",
      "Conflicting Classification" = "#ae017e",
      "Likely Pathogenic" = "#fdae61",
      "Pathogenic" = "#d7191c"
    ),
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE, labels = function(x) str_wrap(x, width = 12)) + 
  theme_minimal(base_size = 14) +
  labs(
    title = "Collated Unique KCNH2 Missense Variant Classification",
    x = "Consensus Classification Status",
    y = "Number of Unique Variants"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.line.x = element_line(color = "black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
  )

print(p2)
ggsave("Figure_02_Overall_Unique_Missense_Classifications.png", plot = p2, width = 10, height = 6, dpi = 300)


## Missense variants with conflicting classifications (Figure 3 Grid Map)
#######################################################################
df_conflicts_base <- df_missense_cleaned %>%
  filter(Classification != "Classification not provided")

conflicting_variants_keys <- df_conflicts_base %>%
  group_by(Final_cnomen) %>%
  summarise(num_unique_classes = n_distinct(Classification), .groups = 'drop') %>%
  filter(num_unique_classes > 1) %>%
  pull(Final_cnomen)

df_conflicts_plot <- df_conflicts_base %>%
  filter(Final_cnomen %in% conflicting_variants_keys) %>%
  mutate(Clean_Variant_Name = gsub(".*:p\\.", "p.", Final_pnomen)) %>%
  group_by(Clean_Variant_Name, Lab) %>%
  summarise(Classification = paste(unique(Classification), collapse = ", "), .groups = 'drop')

class_levels <- c("Likely Benign", "VUS", "Risk Allele", "Likely Pathogenic", "Pathogenic")
df_conflicts_plot$Classification = factor(df_conflicts_plot$Classification, levels = class_levels)

conflict_plot <- ggplot(df_conflicts_plot, aes(x = Lab, y = Clean_Variant_Name, fill = Classification)) +
  geom_tile(color = "white", size = 0.8) +
  scale_fill_manual(
    name = "Laboratory Classification",
    values = c(
      "Likely Benign" = "#abd9e9",
      "VUS" = "#ffffbf",
      "Risk Allele" = "#b2abd2",
      "Likely Pathogenic" = "#fdae61",
      "Pathogenic" = "#d7191c"
    ),
    na.value = "grey95", 
    drop = FALSE
  ) +
  theme_minimal(base_size = 13) +
  labs(
    title = "KCNH2 Missense Variants Without Consensus Classification", 
    x = "Submitting Laboratory",
    y = "Protein Nomenclature"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold", color = "black"),
    axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, family = "mono", color = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

print(conflict_plot)
ggsave("Figure_03_Conflicting_Missense_Variants_Grid_Map.png", plot = conflict_plot, width = 11, height = 7, dpi = 300)


### Variant functional data. 
############################
df_func_breakdown <- Variants_combined_with_APC_functional_data_unique %>%
  filter(
    grepl("^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$", Final_pnomen) &
      !grepl("^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)", Final_pnomen) &
      grepl("^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$", Final_cnomen)
  ) %>%
  mutate(
    Raw_Class = trimws(Lab.variant.classification),
    Clean_Class = case_when(
      grepl("\\|", Raw_Class) ~ "Conflicting Classification",
      Raw_Class %in% c("VUS", "Variant of uncertain significance", "Uncertain Significance") ~ "VUS",
      Raw_Class %in% c("Likely Pathogenic", "Likely pathogenic", "LP") ~ "Likely Pathogenic",
      Raw_Class %in% c("Pathogenic", "pathogenic", "P") ~ "Pathogenic",
      Raw_Class %in% c("Likely Benign", "Likely benign", "LB") ~ "Likely Benign",
      Raw_Class %in% c("Benign", "benign", "B") ~ "Benign",
      TRUE ~ "Conflicting Classification"
    )
  ) %>%
  mutate(
    Has_Match = !is.na(Match_in_APC_functional_df) & (Match_in_APC_functional_df == TRUE | Match_in_APC_functional_df == "True"),
    Has_Z_Score = Final_cnomen %in% Functional_evidence_result$Final_cnomen,
    Functional_Tier = case_when(
      !Has_Match ~ "No APC data",
      Has_Match & !Has_Z_Score ~ "APC data without Z score",
      Has_Match & Has_Z_Score ~ "APC data with Z score"
    )
  ) %>%
  filter(Clean_Class %in% c("Benign", "Likely Benign", "VUS", "Conflicting Classification", "Likely Pathogenic", "Pathogenic"))

df_func_breakdown$Clean_Class <- factor(df_func_breakdown$Clean_Class, levels = c("Benign", "Likely Benign", "VUS", "Conflicting Classification", "Likely Pathogenic", "Pathogenic"))
df_func_breakdown$Functional_Tier <- factor(df_func_breakdown$Functional_Tier, levels = c("No APC data", "APC data without Z score", "APC data with Z score"))

df_matrix_counts <- df_func_breakdown %>%
  group_by(Clean_Class, Functional_Tier) %>%
  summarise(Count = n(), .groups = 'drop')

p4_matrix <- ggplot(df_matrix_counts, aes(x = Functional_Tier, y = Clean_Class, fill = Count)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = Count), fontface = "bold", size = 5, 
            color = ifelse(df_matrix_counts$Count > 35, "white", "black")) +
  scale_fill_gradient(low = "#f5fafd", high = "#08519c", name = "Variant Count") +
  theme_minimal(base_size = 14) +
  labs(
    title = "KCNH2 Missense Variants Functional Data Status",
    x = "Automated Patch-Clamp Curation Pipeline Status",
    y = "Laboratory Classification"
  ) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "right"
  )

print(p4_matrix)
ggsave("Figure_04_Functional_Curation_Matrix.png", plot = p4_matrix, width = 10, height = 5, dpi = 300)


## Z score distribution for VUS variants (Figure 5)
####################################################
df_vus_functional <- Variants_combined_with_APC_functional_data_unique %>%
  filter(
    grepl("^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$", Final_pnomen) &
      !grepl("^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)", Final_pnomen) &
      grepl("^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$", Final_cnomen)
  ) %>%
  filter(grepl("VUS|uncertain", Lab.variant.classification, ignore.case = TRUE)) %>%
  inner_join(Functional_evidence_result, by = c("Final_cnomen", "Final_pnomen")) %>%
  mutate(
    Raw_Class = trimws(Lab.variant.classification),
    Status_Group = if_else(grepl("\\|", Raw_Class), "Conflicting Classification", "Consensus VUS"),
    Assigned_Code = case_when(
      PS3_Strong == "TRUE" | PS3_Strong == TRUE ~ "PS3_Strong",
      PS3_Mod == "TRUE"    | PS3_Mod == TRUE    ~ "PS3_Moderate",
      PS3_Supp == "TRUE"   | PS3_Supp == TRUE   ~ "PS3_Supporting",
      BS3_Supp == "TRUE"   | BS3_Supp == TRUE   ~ "BS3_Supporting",
      BS3_Mod == "TRUE"    | BS3_Mod == TRUE    ~ "BS3_Moderate",
      TRUE ~ "No Criteria Applied"
    ),
    Short_p_nomen = gsub("^.*:p\\.", "p.", Final_pnomen)
  )

df_vus_functional$Short_p_nomen <- reorder(df_vus_functional$Short_p_nomen, -df_vus_functional$Z_score_used)
df_vus_functional$Assigned_Code <- factor(df_vus_functional$Assigned_Code, levels = c("BS3_Moderate", "BS3_Supporting", "PS3_Supporting", "PS3_Moderate", "PS3_Strong"))
df_vus_functional$Status_Group <- factor(df_vus_functional$Status_Group, levels = c("Consensus VUS", "Conflicting Classification"))

p5_vus_distribution <- ggplot(df_vus_functional, aes(x = Short_p_nomen, y = Z_score_used, color = Assigned_Code, shape = Status_Group)) +
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = -Inf, ymax = -4, fill = "#d7191c", alpha = 0.05) + 
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = -4, ymax = -3, fill = "#fdae61", alpha = 0.05) +   
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = -3, ymax = -2, fill = "#fee090", alpha = 0.06) +   
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = -2, ymax = -1, fill = "#abd9e9", alpha = 0.04) +   
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = -1, ymax = 1,  fill = "#2c7bb6", alpha = 0.03) +   
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 1,  ymax = 2,  fill = "#abd9e9", alpha = 0.04) +   
  geom_hline(yintercept = c(-4, -3, -2, -1, 1, 2), linetype = "dashed", color = "grey55", size = 0.5) +
  geom_point(size = 3.5, alpha = 0.95) +
  scale_color_manual(
    name = "Assigned ACGS Code",
    values = c("BS3_Moderate" = "#2c7bb6", "BS3_Supporting" = "#abd9e9", "PS3_Supporting" = "#e0b034", "PS3_Moderate" = "#fdae61", "PS3_Strong" = "#d7191c")
  ) +
  scale_shape_manual(
    name = "Laboratory Classification Status",
    values = c("Consensus VUS" = 16, "Conflicting Classification" = 17)
  ) +
  scale_y_continuous(limits = c(-8, 3), breaks = seq(-8, 3, by = 2)) +
  theme_minimal(base_size = 14) +
  labs(title = "Functional Evidence Strength for KCNH2 VUS Missense Variants", x = "Variant Protein Nomenclature", y = "Z-score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black"),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15))
  )

print(p5_vus_distribution)
ggsave("Figure_05_VUS_Functional_Z_Score_Distribution.png", plot = p5_vus_distribution, width = 12, height = 7, dpi = 300)


### Automated criteria for variants: Ingest & Keys
###################################################
df_vus_core <- df_vus_functional %>% mutate(short_c = gsub("^NM_000238(\\.[0-9]+)?:", "", Final_cnomen))
df_gnomad_clean <- Allele_freq_criterion_result %>% mutate(short_c = gsub("^NM_000238(\\.[0-9]+)?:", "", Match_key)) %>% distinct(short_c, .keep_all = TRUE)
df_revel_clean <- REVEL_evidence_result %>% mutate(short_c = gsub("^NM_000238(\\.[0-9]+)?:", "", Variant_cnomen)) %>% distinct(short_c, .keep_all = TRUE)
df_splice_clean <- SpliceAI_evidence_result %>% mutate(short_c = gsub("^NM_000238(\\.[0-9]+)?:", "", Variant_cnomen)) %>% distinct(short_c, .keep_all = TRUE)
df_pm1_clean <- PM1_domain_criterion_results %>% mutate(short_c = gsub("^NM_000238(\\.[0-9]+)?:", "", Final_cnomen)) %>% distinct(short_c, .keep_all = TRUE)

df_automated_dashboard <- df_vus_core %>%
  select(short_c, Short_p_nomen, Z_score_used, Status_Group) %>%
  left_join(df_gnomad_clean %>% select(short_c, FAF_FAFpop, LowerCI_FAFpop, UpperCI_FAFpop, BA1, BS1, PM2_Supporting), by = "short_c") %>%
  left_join(df_revel_clean  %>% select(short_c, REVEL_score), by = "short_c") %>%
  left_join(df_splice_clean %>% select(short_c, Max_Delta_Score), by = "short_c")

df_automated_dashboard$REVEL_score <- as.numeric(df_automated_dashboard$REVEL_score)
df_automated_dashboard$Short_p_nomen <- reorder(df_automated_dashboard$Short_p_nomen, -df_automated_dashboard$Z_score_used)

df_automated_dashboard <- df_automated_dashboard %>%
  mutate(
    FAF_FAFpop     = if_else(is.na(FAF_FAFpop), 1e-7, FAF_FAFpop),
    LowerCI_FAFpop = if_else(is.na(LowerCI_FAFpop), 1e-7, LowerCI_FAFpop),
    UpperCI_FAFpop = if_else(is.na(UpperCI_FAFpop), 1e-7, UpperCI_FAFpop),
    Pop_Criterion = case_when(
      BA1 == "TRUE" | BA1 == TRUE ~ "BA1",
      BS1 == "TRUE" | BS1 == TRUE ~ "BS1",
      PM2_Supporting == "TRUE" | PM2_Supporting == TRUE ~ "PM2",
      TRUE ~ "None (Intermediate Freq)"
    ),
    Pop_Criterion = factor(Pop_Criterion, levels = c("BA1", "BS1", "None (Intermediate Freq)", "PM2"))
  )


## VUS variants GnomAD Data (Figure 6)
######################################
p6_gnomad_standalone <- ggplot(df_automated_dashboard, aes(x = Short_p_nomen, y = FAF_FAFpop, color = Pop_Criterion)) +
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0.001, ymax = Inf, fill = "#08519c", alpha = 0.05) + 
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0.00005, ymax = 0.001, fill = "#2c7bb6", alpha = 0.05) + 
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0, ymax = 0.00001, fill = "#d7191c", alpha = 0.05) +   
  geom_hline(yintercept = c(0.001, 0.00005, 0.00001), linetype = "dashed", color = "grey30", size = 0.6) +
  geom_errorbar(aes(ymin = pmax(LowerCI_FAFpop, 1e-7), ymax = UpperCI_FAFpop), width = 0.2, color = "grey45") +
  geom_point(aes(shape = Status_Group), size = 3.5, alpha = 0.95, stroke = 0.8) +
  scale_color_manual(
    name = "Triggered Population Criterion",
    values = c("BA1" = "#08519c", "BS1" = "#2c7bb6", "None (Intermediate Freq)" = "#969696", "PM2" = "#d7191c"),
    drop = FALSE
  ) +
  scale_shape_manual(name = "Laboratory Status", values = c("Consensus VUS" = 16, "Conflicting Classification" = 17)) +
  scale_x_discrete() +
  scale_y_log10(limits = c(5e-8, 1e-2), 
                breaks = c(1e-7, 1e-6, 1e-5, 5e-5, 1e-3, 1e-2),
                labels = c("0 / Not detected", "0.000001", "0.00001", "0.00005", "0.001", "0.01")) +
  theme_minimal(base_size = 14) +
  labs(title = "GnomAD Filtered Allele Frequency Spectrum for KCNH2 VUS Cohort", x = "Variant Protein Nomenclature", y = "GnomAD Filtered Allele Frequency") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black"),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15))
  )

print(p6_gnomad_standalone)
ggsave("Figure_06_GnomAD_Population_Frequency_Spectrum.png", plot = p6_gnomad_standalone, width = 12, height = 7, dpi = 300)


## VUS variants REVEL Data (Figure 7)
#####################################
df_revel_standalone <- df_vus_core %>%
  select(short_c, Short_p_nomen, Z_score_used, Status_Group) %>%
  left_join(df_revel_clean %>% select(short_c, REVEL_score), by = "short_c") %>%
  mutate(
    REVEL_score = as.numeric(REVEL_score),
    REVEL_criterion = case_when(
      REVEL_score >= 0.70 ~ "PP3",
      REVEL_score <= 0.40 ~ "BP4",
      TRUE ~ "None (Intermediate Score)"
    ),
    REVEL_criterion = factor(REVEL_criterion, levels = c("BP4", "None (Intermediate Score)", "PP3"))
  )

df_revel_standalone$Short_p_nomen <- reorder(df_revel_standalone$Short_p_nomen, -df_revel_standalone$Z_score_used)

p7_revel_standalone <- ggplot(df_revel_standalone, aes(x = Short_p_nomen, y = REVEL_score, color = REVEL_criterion)) +
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0.70, ymax = 1.00, fill = "#d7191c", alpha = 0.05) + 
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0.00, ymax = 0.40, fill = "#2c7bb6", alpha = 0.05) + 
  geom_hline(yintercept = c(0.40, 0.70), linetype = "dashed", color = "grey30", size = 0.6) +
  geom_point(aes(shape = Status_Group), size = 3.5, alpha = 0.95, stroke = 0.8) +
  scale_color_manual(
    name = "Triggered REVEL Criterion",
    values = c("BP4" = "#2c7bb6", "None (Intermediate Score)" = "#969696", "PP3" = "#d7191c"),
    drop = FALSE
  ) +
  scale_shape_manual(name = "Laboratory Status", values = c("Consensus VUS" = 16, "Conflicting Classification" = 17)) +
  scale_x_discrete() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), labels = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0")) +
  theme_minimal(base_size = 14) +
  labs(title = "REVEL Ensemble Pathogenicity Spectrum for KCNH2 VUS Cohort", x = "Variant Protein Nomenclature", y = "REVEL Score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black"),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15))
  )

print(p7_revel_standalone)
ggsave("Figure_07_REVEL_Pathogenicity_Scores.png", plot = p7_revel_standalone, width = 12, height = 7, dpi = 300)


## VUS variants SpliceAI Data (Figure 8)
########################################
df_splice_standalone <- df_vus_core %>%
  select(short_c, Short_p_nomen, Z_score_used, Status_Group) %>%
  left_join(df_splice_clean %>% select(short_c, Max_Delta_Score), by = "short_c") %>%
  mutate(Max_Delta_Score = as.numeric(Max_Delta_Score))

df_splice_standalone$Short_p_nomen <- reorder(df_splice_standalone$Short_p_nomen, -df_splice_standalone$Z_score_used)

p8_spliceai_standalone <- ggplot(df_splice_standalone, aes(x = Short_p_nomen, y = Max_Delta_Score)) +
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0.20, ymax = 0.30, fill = "#d7191c", alpha = 0.05) + 
  annotate("rect", xmin = 0.4, xmax = 38.6, ymin = 0.00, ymax = 0.10, fill = "#2c7bb6", alpha = 0.05) + 
  geom_hline(yintercept = c(0.10, 0.20), linetype = "dashed", color = "grey30", size = 0.6) +
  annotate("text", x = 38.5, y = 0.20, label = "PP3 (> 0.2)", fontface = "bold", size = 4, color = "#d7191c", hjust = 1, vjust = -0.4) +
  annotate("text", x = 38.5, y = 0.10, label = "BP4 (< 0.1)", fontface = "bold", size = 4, color = "#2c7bb6", hjust = 1, vjust = -0.4) +
  geom_point(aes(shape = Status_Group), size = 3.5, alpha = 0.95, stroke = 0.8, color = "black") +
  scale_shape_manual(name = "Laboratory Status", values = c("Consensus VUS" = 16, "Conflicting Classification" = 17)) +
  scale_x_discrete() +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 0.3, 0.1), labels = c("0.0", "0.1", "0.2", "0.3")) +
  theme_minimal(base_size = 14) +
  labs(title = "SpliceAI Splicing Disturbance Spectrum for KCNH2 VUS Cohort", x = "Variant Protein Nomenclature", y = "Maximum Delta Score") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black"),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15))
  )

print(p8_spliceai_standalone)
ggsave("Figure_08_SpliceAI_Splicing_Disturbance_Scores.png", plot = p8_spliceai_standalone, width = 12, height = 7, dpi = 300)


### Unified Computational Predictor Dashboard Panel (Figure 9)
###############################################################
p6_panel <- p6_gnomad_standalone +
  labs(title = "A) GnomAD Filtered Allele Frequency Spectrum", x = NULL) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )

p7_panel <- p7_revel_standalone +
  labs(title = "B) REVEL Pathogenicity Spectrum", x = NULL) +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )

p8_panel <- p8_spliceai_standalone +
  labs(title = "C) SpliceAI Splicing Disturbance Spectrum", x = "Variant Protein Nomenclature") +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "bold", color = "black", size = 8.5),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.position = "right"
  )

composite_portrait_dashboard <- p6_panel / p7_panel / p8_panel + plot_layout(heights = c(1, 1, 1))

print(composite_portrait_dashboard)
ggsave("Figure_09_Unified_Computational_Predictor_Dashboard.png", 
       plot = composite_portrait_dashboard, width = 9.5, height = 12.5, dpi = 300)


## VUS variants PM1 Structural Mapping (Figure 10)
##################################################
if(!require(ggrepel)) install.packages("ggrepel")
library(ggrepel)

df_pm1_map <- df_vus_core %>%
  left_join(df_pm1_clean %>% select(short_c, PM1_mod, PM1_sup, Domain_Comment), by = "short_c") %>%
  mutate(
    Original_Codon = as.numeric(stringr::str_extract(Final_pnomen, "(?<=p\\.[A-Za-z]{3})[0-9]+|(?<=p\\.[A-Z])[0-9]+")),
    Projected_Codon = Original_Codon * 128,
    Structural_Tier = case_when(
      PM1_mod == "TRUE" | PM1_mod == TRUE ~ "PM1_Moderate Hotspot",
      PM1_sup == "TRUE" | PM1_sup == TRUE ~ "PM1_Supporting Hotspot",
      TRUE ~ "Outside Hotspot Domains"
    )
  ) %>%
  # FIXED: Swapped problematic negative operator with logical logical NOT expression
  filter(!is.na(Projected_Codon))

exon_coordinates <- data.frame(
  Exon  = paste0("E", 1:15),
  Start = c(1,   26,  103, 158, 306, 377, 520, 649, 716, 800, 865, 898, 989,  1051, 1111) * 128,
  End   = c(26,  103, 158, 306, 376, 519, 649, 715, 800, 864, 898, 989, 1051, 1110, 1159) * 128
)

p9_pm1_map <- ggplot() +
  annotate("rect", xmin = 1*128,   xmax = 130*128, ymin = 0.00, ymax = 0.06, fill = "#d7191c", alpha = 0.20) + 
  annotate("rect", xmin = 404*128, xmax = 659*128, ymin = 0.00, ymax = 0.06, fill = "#d7191c", alpha = 0.20) + 
  annotate("rect", xmin = 742*128, xmax = 842*128, ymin = 0.00, ymax = 0.06, fill = "#fdae61", alpha = 0.20) + 
  annotate("text", x = 65*128,  y = 0.03, label = "N-terminus Cluster (Mod)", fontface = "bold", size = 2.4, color = "#b30000") +
  annotate("text", x = 531*128, y = 0.03, label = "Transmembrane / Linker / Pore (Mod)", fontface = "bold", size = 2.5, color = "#b30000") +
  annotate("text", x = 792*128, y = 0.03, label = "C-terminus cNBD (Sup)", fontface = "bold", size = 2.4, color = "#d95f02") +
  annotate("text", x = 1*128,    y = -0.04, label = "1",   fontface = "bold", size = 3, color = "black", hjust = 0) +
  annotate("text", x = 130*128,  y = -0.04, label = "130", fontface = "bold", size = 3, color = "black", hjust = 1) +
  annotate("text", x = 404*128,  y = -0.04, label = "404", fontface = "bold", size = 3, color = "black", hjust = 0) +
  annotate("text", x = 659*128,  y = -0.04, label = "659", fontface = "bold", size = 3, color = "black", hjust = 1) +
  annotate("text", x = 742*128,  y = -0.04, label = "742", fontface = "bold", size = 3, color = "black", hjust = 0) +
  annotate("text", x = 842*128,  y = -0.04, label = "842", fontface = "bold", size = 3, color = "black", hjust = 1) +
  geom_rect(data = exon_coordinates, aes(xmin = Start, xmax = End, ymin = 0.060, ymax = 0.135), fill = "grey93", color = "black", size = 0.4) +
  geom_text(data = exon_coordinates, aes(x = (Start + End)/2, y = 0.0975, label = Exon), fontface = "bold", size = 2.4, color = "black") +
  geom_segment(data = df_pm1_map, aes(x = Projected_Codon, xend = Projected_Codon, y = 0.136, yend = 0.180), color = "grey50", size = 0.4, linetype = "solid") +
  geom_point(data = df_pm1_map, aes(x = Projected_Codon, y = 0.180, color = Structural_Tier, shape = Status_Group), size = 3.5, stroke = 0.8) +
  geom_text_repel(
    data = df_pm1_map, aes(x = Projected_Codon, y = 0.195, label = Short_p_nomen),
    angle = 90, hjust = 0, vjust = 0.5, fontface = "bold", size = 3.2, color = "black",
    direction = "x", nudge_y = 0.015, box.padding = 0.12, point.padding = 0.06, 
    segment.color = "grey60", segment.size = 0.35, max.overlaps = Inf      
  ) +
  scale_color_manual(name = "Triggered PM1 Criterion:", values = c("PM1_Moderate Hotspot" = "#d7191c", "PM1_Supporting Hotspot" = "#fdae61", "Outside Hotspot Domains" = "#969696")) +
  scale_shape_manual(name = "Laboratory Classification:", values = c("Consensus VUS" = 16, "Conflicting Classification" = 17)) +
  scale_x_continuous(limits = c(1 * 128, 1160 * 128), expand = c(0.015, 0.015)) +
  scale_y_continuous(limits = c(-0.08, 0.58)) +
  theme_minimal(base_size = 14) +
  labs(title = "KCNH2 Linear Structural Mapping & PM1 Functional Hotspot Allocation", x = "KCNH2 Amino Acid Codon Coordinate Map (Codons 1 to 1159)", y = NULL) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 12, b = 10)),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 15)),
    legend.position = "bottom", legend.box = "horizontal", legend.title = element_text(face = "bold", size = 11), legend.margin = margin(t = 10)
  ) +
  guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))

print(p9_pm1_map)
ggsave("Figure_10_KCNH2_PM1_Structural_Map.png", plot = p9_pm1_map, width = 20, height = 5.6, dpi = 300)


### Master 38 VUS Variant Summary Table for Thesis Layout
#########################################################
df_matrix <- read.csv("Combined_criteria_matrix_results - EDIT.csv", stringsAsFactors = FALSE)

criteria_cols <- c("PVS1", "PS1", "PS2", "PS3", "PS4", "PM1", "PM4", "PM5", "PM6", 
                   "PP1", "PP2", "PP3", "PM2", "BA1", "BS1", "BS3", "BS4", "BP2", "BP4", "BP7")

format_applied_criteria <- function(df_row) {
  applied <- c()
  for (col in criteria_cols) {
    val <- df_row[[col]]
    if (!is.na(val) && val != "" && val != 0) {
      strength <- switch(trimws(as.character(val)),
                         "VStr" = "very_strong", "Str"  = "strong", "Mod"  = "moderate", "Sup"  = "supporting",
                         tolower(trimws(as.character(val))))
      applied <- c(applied, paste0(col, "_", strength))
    }
  }
  if (length(applied) == 0) return("None")
  return(paste(applied, collapse = ", "))
}

master_vus_table <- df_matrix %>%
  filter(grepl("VUS", Collated_classification, ignore.case = TRUE)) %>%
  rowwise() %>%
  mutate(
    Criteria_Applied = format_applied_criteria(pick(all_of(criteria_cols))),
    Clean_p_nomen = gsub(".*:p\\.", "p.", Final_pnomen)
  ) %>%
  ungroup() %>%
  select(
    `cDNA Nomenclature`    = Variant_cnomen,
    `Protein Nomenclature`  = Clean_p_nomen,
    `ACGS Criteria Applied` = Criteria_Applied,
    `Total Points`           = Total_ACGS_Points,
    `Reclassification`      = Final_classification 
  ) %>%
  arrange(as.numeric(str_extract(`cDNA Nomenclature`, "[0-9]+")))

write.csv(master_vus_table, "Table_1_VUS_Variant_Summary.csv", row.names = FALSE)


### Live-Calculated Gradient Alluvial Curation Flow Charts (ACGS 2024 Compliant)
###############################################################################
# 1. Custom Sigmoid Gradient Mesh Generator Function
generate_gradient_ribbons <- function(data, col_left, col_right, x_left, x_right, order_left, order_right, colors_left, colors_right) {
  flows <- data %>%
    group_by(Left = .data[[col_left]], Right = .data[[col_right]]) %>%
    summarise(Freq = n(), .groups = 'drop') %>%
    filter(Freq > 0)
  
  flows <- flows %>% arrange(match(Left, order_left), match(Right, order_right))
  flows$y_left_min <- c(0, cumsum(flows$Freq)[-nrow(flows)])
  flows$y_left_max <- cumsum(flows$Freq)
  
  flows <- flows %>% arrange(match(Right, order_right), match(Left, order_left))
  flows$y_right_min <- c(0, cumsum(flows$Freq)[-nrow(flows)])
  flows$y_right_max <- cumsum(flows$Freq)
  
  steps <- 100
  sig <- 1 / (1 + exp(-seq(-6, 6, length.out = steps)))
  sig <- (sig - min(sig)) / (max(sig) - min(sig)) 
  
  slices <- list()
  for (i in 1:nrow(flows)) {
    row <- flows[i, ]
    hex_left  <- colors_left[as.character(row$Left)]
    hex_right <- colors_right[as.character(row$Right)]
    flow_colors <- colorRampPalette(c(hex_left, hex_right))(steps)
    
    for (j in 1:(steps - 1)) {
      t1 <- sig[j]; t2 <- sig[j + 1]
      xl <- x_left + (j - 1) / (steps - 1) * (x_right - x_left)
      xr <- x_left + j / (steps - 1) * (x_right - x_left)
      yl1 <- row$y_left_min + (row$y_right_min - row$y_left_min) * t1
      yl2 <- row$y_left_max + (row$y_right_max - row$y_left_max) * t1
      yr1 <- row$y_left_min + (row$y_right_min - row$y_left_min) * t2
      yr2 <- row$y_left_max + (row$y_right_max - row$y_left_max) * t2
      
      slices[[length(slices) + 1]] <- data.frame(
        flow_id = paste0(i, "_", j), x = c(xl, xr, xr, xl),
        y = c(yl1, yr1, yr2, yl2), fill = flow_colors[j], stringsAsFactors = FALSE
      )
    }
  }
  return(do.call(rbind, slices))
}

# 2. Custom Stratum Block Vector Calculator Function
generate_stratum_blocks <- function(data, column, x_pos, node_order, width = 0.35) {
  total_variants <- nrow(data)
  counts <- data %>%
    group_by(Node = .data[[column]]) %>%
    summarise(Freq = n(), .groups = 'drop') %>%
    filter(Freq > 0) %>%
    arrange(match(Node, node_order)) %>%
    mutate(
      Pct = (Freq / total_variants) * 100,
      Label = case_when(
        Freq == 1 ~ paste0(Node, " (n=", Freq, ", ", sprintf("%.1f", Pct), "%)"),
        TRUE      ~ paste0(str_wrap(Node, width = 12), "\n(n=", Freq, ", ", sprintf("%.1f", Pct), "%)")
      ),
      y_pos = (cumsum(Freq) + c(0, cumsum(Freq)[-n()])) / 2,
      f_size = case_when(Freq == 1 ~ 2.0, TRUE ~ 3.2)
    )
  counts$ymax <- cumsum(counts$Freq)
  counts$ymin <- c(0, counts$ymax[-nrow(counts)])
  counts$xmin <- x_pos - width/2
  counts$xmax <- x_pos + width/2
  return(counts)
}

# 3. Sorting Dimensions and Color Palettes
order_before  <- c("Consensus VUS", "Conflicting Classification")
order_after   <- c("Benign", "Likely Benign", "VUS", "Likely Pathogenic", "Pathogenic")
colors_before <- c("Consensus VUS" = "#4daf4a", "Conflicting Classification" = "#984ea3")
colors_after  <- c("Benign" = "#2c7bb6", "Likely Benign" = "#abd9e9", "VUS" = "#ffffbf", "Likely Pathogenic" = "#fdae61", "Pathogenic" = "#d7191c")

### Execute Direct Matrix Calculations (Strict ACGS 2024 Framework)
##################################################################
df_matrix_calculated <- df_matrix %>%
  filter(grepl("VUS", Collated_classification, ignore.case = TRUE)) %>%
  mutate(
    Before = case_when(
      grepl("\\|", Collated_classification) | grepl("conflict", Collated_classification, ignore.case = TRUE) ~ "Conflicting Classification",
      TRUE ~ "Consensus VUS"
    ),
    PS3_pts = case_when(trimws(as.character(PS3)) == "Str" ~ 4, trimws(as.character(PS3)) == "Mod" ~ 2, trimws(as.character(PS3)) == "Sup" ~ 1, TRUE ~ 0),
    BS3_pts = case_when(trimws(as.character(BS3)) == "Str" ~ -4, trimws(as.character(BS3)) == "Mod" ~ -2, trimws(as.character(BS3)) == "Sup" ~ -1, TRUE ~ 0),
    PS4_pts = case_when(trimws(as.character(PS4)) == "Str" ~ 4, trimws(as.character(PS4)) == "Mod" ~ 2, trimws(as.character(PS4)) == "Sup" ~ 1, TRUE ~ 0),
    PP1_pts = case_when(trimws(as.character(PP1)) == "Str" ~ 4, trimws(as.character(PP1)) == "Mod" ~ 2, trimws(as.character(PP1)) == "Sup" ~ 1, TRUE ~ 0),
    
    Criteria_Count = (PVS1!=0) + (PS1!=0) + (PS2!=0) + (PS3!=0) + (PS4!=0) + (PM1!=0) + (PM4!=0) + (PM5!=0) + (PM6!=0) + (PP1!=0) + (PP2!=0) + (PP3!=0) +
      (BA1!=0) + (BS1!=0) + (BS3!=0) + (BS4!=0) + (BP2!=0) + (BP4!=0) + (BP7!=0),
    
    Classification_Revised = case_when(
      (Criteria_Count < 2 & BA1 != "VStr") ~ "VUS",
      Total_ACGS_Points >= 10 ~ "Pathogenic", Total_ACGS_Points >= 6 ~ "Likely Pathogenic",
      Total_ACGS_Points >= 0 ~ "VUS", Total_ACGS_Points >= -5 ~ "Likely Benign",
      Total_ACGS_Points <= -6 ~ "Benign", TRUE ~ "VUS"
    ),
    Net_Points_No_Functional = Total_ACGS_Points - PS3_pts - BS3_pts,
    Classification_No_Functional = case_when(
      (Criteria_Count < 2 & BA1 != "VStr") ~ "VUS",
      Net_Points_No_Functional >= 10 ~ "Pathogenic", Net_Points_No_Functional >= 6 ~ "Likely Pathogenic",
      Net_Points_No_Functional >= 0 ~ "VUS", Net_Points_No_Functional >= -5 ~ "Likely Benign",
      Net_Points_No_Functional <= -6 ~ "Benign", TRUE ~ "VUS"
    ),
    Net_Points_No_Pedigree = Total_ACGS_Points - PS4_pts - PP1_pts,
    Classification_No_Pedigree = case_when(
      (Criteria_Count < 2 & BA1 != "VStr") ~ "VUS",
      Net_Points_No_Pedigree >= 10 ~ "Pathogenic", Net_Points_No_Pedigree >= 6 ~ "Likely Pathogenic",
      Net_Points_No_Pedigree >= 0 ~ "VUS", Net_Points_No_Pedigree >= -5 ~ "Likely Benign",
      Net_Points_No_Pedigree <= -6 ~ "Benign", TRUE ~ "VUS"
    )
  )

## Plot Rendering Pipeline (Figure 13 Composite)
################################################
generate_sankey_object <- function(data, target_col, title_txt, right_label, add_triangle = FALSE) {
  mesh <- generate_gradient_ribbons(data, "Before", target_col, 1, 2, order_before, order_after, colors_before, colors_after)
  b1   <- generate_stratum_blocks(data, "Before", 1, order_before)
  b2   <- generate_stratum_blocks(data, target_col, 2, order_after)
  
  p <- ggplot() +
    geom_polygon(data = mesh, aes(x = x, y = y, group = flow_id, fill = fill), alpha = 0.60) +
    geom_rect(data = b1, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey98", color = "black", size = 0.4) +
    geom_rect(data = b2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey98", color = "black", size = 0.4) +
    
    geom_text(data = b1, aes(x = 1, y = y_pos, label = Label, size = f_size), fontface = "bold", lineheight = 0.95, vjust = 0.5) +
    geom_text(data = b2, aes(x = 2, y = y_pos, label = Label, size = f_size), fontface = "bold", lineheight = 0.95, vjust = 0.5) +
    
    annotate("text", x = 1, y = 39.7, label = "Submitted Lab Classification", fontface = "bold", size = 3.8, color = "black") +
    annotate("text", x = 2, y = 39.7, label = right_label, fontface = "bold", size = 3.8, color = "black")
  
  if (add_triangle) {
    benign_y <- b2$y_pos[b2$Node == "Benign"]
    if (length(benign_y) > 0) {
      p <- p + annotate("point", x = 2.22, y = benign_y, shape = 17, size = 3.5, color = "black")
    }
  }
  
  p <- p +
    scale_fill_identity() +
    scale_size_identity() + 
    scale_y_reverse(limits = c(40.5, -1.5), expand = c(0,0)) + 
    scale_x_continuous(limits = c(0.65, 2.35), expand = c(0, 0)) +
    theme_minimal(base_size = 12) +
    labs(title = title_txt, x = NULL, y = NULL) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5, margin = margin(t = 10, b = 10))
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# Generate individual plots
fig11 <- generate_sankey_object(df_matrix_calculated, "Classification_Revised", "Reclassification of 38 KCNH2 Missense VUS Variants", "Reclassification", add_triangle = FALSE)
fig12 <- generate_sankey_object(df_matrix_calculated, "Classification_No_Functional", "Variant Reclassification in the Absence of PS3/BS3 Evidence", "Reclassification (No PS3/BS3 Criteria Applied)", add_triangle = TRUE)
fig13 <- generate_sankey_object(df_matrix_calculated, "Classification_No_Pedigree", "Variant Reclassification in the Absence of PS4 or PP1 Evidence", "Reclassification (No PS4 or PP1 Criteria Applied)", add_triangle = FALSE)

# Compile vertical stack
thesis_composite_sankey <- (fig11 / fig12 / fig13) + 
  plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A", tag_suffix = ")") & 
  theme(
    plot.margin = margin(t = 12, b = 12, l = 5, r = 25), 
    plot.tag = element_text(face = "bold", size = 16)
  )

# Output image
print(thesis_composite_sankey)
ggsave("Figure_13_Reclassification_Composite_Sankey.png", plot = thesis_composite_sankey, width = 11, height = 16.5, dpi = 300)
