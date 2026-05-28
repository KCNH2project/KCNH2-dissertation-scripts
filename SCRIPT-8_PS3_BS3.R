### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load in required libraries.
###############################
library(lubridate)
library(dplyr)
library(stringr)

### Load in variant data.
#########################
Collated_variant <- read.csv("Results_variant_recoder_protein_df_final.csv")

### Load in the updated Oxford classifications from script 5. 
#############################################################
Oxford_classifications <- read.csv("Oxford_classifications_for_all_lab_variants - EDIT.csv")

### Load in functional data. 
############################
Merged_functional_df <- read.csv("Merged_functional_data.csv")

### Edit the collated variants to take into account the updated Oxford classifications.
#######################################################################################
## Create a temporary clean nomenclature column for precise cross-referencing
Collated_variant$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", Collated_variant$Final_cnomen)

## Left join the relevant classification and date columns from the Oxford file
Collated_variant <- Collated_variant %>%
  left_join(Oxford_classifications %>% select(Variant_cnomen, Final_classification, Classified_date), 
            by = c("Match_cnomen" = "Variant_cnomen"))

## Standardize the classification terminology layout style to match other labs
Collated_variant <- Collated_variant %>%
  mutate(
    Final_classification = case_when(
      Final_classification == "P"   ~ "Pathogenic",
      Final_classification == "LP"  ~ "Likely pathogenic",
      Final_classification == "LB"  ~ "Likely benign",
      Final_classification == "B"   ~ "Benign",
      TRUE                          ~ Final_classification  # Retains VUS, Risk allele, or unchanged metrics
    )
  )

## Identify rows belonging to the Oxford lab where a valid classification exists in the Oxford file
Oxford_rows <- Collated_variant$Lab == "Oxford" & 
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

### Remove rows with empty classifications. 
###########################################
## This was the case for variant's provided by Manchester, unable to use for further analysis due to lack of
## classifications. 
Collated_variant <- Collated_variant %>%
  filter(!(Lab == "Manchester" & (is.na(Lab.variant.classification) | trimws(Lab.variant.classification) == "")))

### Edit collated dataset to account for duplicates in Brompton data. 
#####################################################################
## Brompton has included multiple classifications for some variants. In such instances, select only the variant with the
## most recent classification date. If one/both/neither have classification dates, retain all, as cannot be sure which
## classification was most recent. If all the remaining duplicates have the same classification and classification
## date, keep only one row. If the remaining duplicates have the same classification date but differing classifications,
## keep the entry with the highest lab number (i.e. the most recently tested and classified patient).

# Ensure the temporary matching column exists for grouping
Collated_variant$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", Collated_variant$Final_cnomen)

# 1. Separate non-Brompton data to ensure it remains completely unaffected
non_brompton_df <- Collated_variant %>% filter(Lab != "Brompton")

# 2. Isolate Brompton data and parse the date and lab number fields for numeric comparison
brompton_df <- Collated_variant %>% 
  filter(Lab == "Brompton") %>%
  mutate(
    # Try parsing common formats. Adjust orders if your dates are formatted differently (e.g., "dmy" or "ymd")
    Parsed_Date = parse_date_time(Classification.date, orders = c("dmy", "ymd", "mdy")),
    # Extract only the numerical digits from the Lab.number column for maximum-value identification
    Numeric_Lab_Num = as.numeric(gsub("[^0-9]", "", Lab.number))
  ) %>%
  mutate(
    # Safely convert missing or unparseable lab numbers to 0 to prevent downstream evaluation issues
    Numeric_Lab_Num = if_else(is.na(Numeric_Lab_Num), 0, Numeric_Lab_Num)
  )

# 3. Clean Brompton data by sequentially evaluating dates, classifications, and lab numbers
brompton_cleaned <- brompton_df %>%
  group_by(Match_cnomen) %>%
  mutate(
    # Determine if any record in the duplicate variant group is missing a classification date
    Any_Date_Missing = any(is.na(Classification.date) | trimws(Classification.date) == ""),
    # Determine the maximum date within the duplicate variant group (ignoring NAs)
    Max_Group_Date = if(all(is.na(Parsed_Date))) as.POSIXct(NA) else max(Parsed_Date, na.rm = TRUE)
  ) %>%
  filter(
    # Step 1: Keep only rows matching the most recent date, OR retain all rows if any date in the group is missing
    Any_Date_Missing | (!is.na(Parsed_Date) & Parsed_Date == Max_Group_Date)
  ) %>%
  mutate(
    # Identify the highest lab number remaining in the filtered variant sub-group
    Max_Lab_Num = max(Numeric_Lab_Num, na.rm = TRUE)
  ) %>%
  filter(
    # Step 2: Keep only the entry (or entries) matching the highest lab number
    Numeric_Lab_Num == Max_Lab_Num
  ) %>%
  # Step 3: If there's an exact tie for date and lab number, keep only the first row to ensure a single unique variant entry
  slice(1) %>%
  ungroup() %>%
  # Remove the helper columns generated for the sequential evaluation calculations
  select(-Parsed_Date, -Numeric_Lab_Num, -Any_Date_Missing, -Max_Group_Date, -Max_Lab_Num)

# 4. Recombine the non-Brompton rows with the cleaned Brompton dataset
Collated_variant <- rbind(non_brompton_df, brompton_cleaned)

# Clean up the temporary matching column
Collated_variant$Match_cnomen <- NULL

## Save updated variant list. 
write.csv(Collated_variant, "Results_variant_recoder_protein_df_updated.csv", row.names = TRUE, quote = TRUE)

### Check if APC data exists for list of lab variants. 
######################################################
## Add protein HGVS to functional protein column. 
Merged_functional_df$`Protein.Change` <- paste0("NP_000229.1:", Merged_functional_df$`Protein.Change`)

## Create a column which determines if there is a match in both data frames, true or false. 
Collated_variant_query_APC_functional_match <- Collated_variant %>%
  mutate(Match_in_APC_functional_df = Final_pnomen %in% Merged_functional_df$`Protein.Change`)

## Merge the data frames using the protein nomenclature. 
Variants_combined_with_APC_functional_data <- Collated_variant_query_APC_functional_match %>%
  left_join(Merged_functional_df, by = c("Final_pnomen" = "Protein.Change"))

## Save and manually review. 
write.csv(Variants_combined_with_APC_functional_data, "Variants_combined_with_APC_functional_data.csv", row.names = TRUE, quote = TRUE)
# Looks correct. 

## Make all variant entries unique.
Variants_combined_with_APC_functional_data_unique <- Variants_combined_with_APC_functional_data %>%
  group_by(Final_cnomen) %>%
  summarise(
    across(
      everything(),
      ~ {
        vals <- trimws(na.omit(.x))              # remove NA + strip whitespace
        vals <- vals[vals != ""]                 # drop empty strings
        if (length(vals) == 0) NA_character_ else paste(unique(vals), collapse = " | ")
      }
    ),
    .groups = "drop"
  )

## Save and manually review. 
write.csv(Variants_combined_with_APC_functional_data_unique, "Variants_combined_with_APC_functional_data_unique.csv", row.names = TRUE, quote = TRUE)

## Pull out variants with no functional data. 
Variants_no_APC_functional_unique <- Variants_combined_with_APC_functional_data_unique[Variants_combined_with_APC_functional_data_unique$Match_in_APC_functional_df == FALSE,]
Variants_with_APC_functional_unique <- Variants_combined_with_APC_functional_data_unique[Variants_combined_with_APC_functional_data_unique$Match_in_APC_functional_df == TRUE,]

### Identify variants with/without Z score.
###########################################

Variants_combined_with_APC_functional_data_unique <- read.csv("Variants_combined_with_APC_functional_data_unique.csv")

## Keep all variants where APC functional data matches, completely irrespective of clinical classification
all_with_apc_data <- Variants_combined_with_APC_functional_data_unique %>%
  filter(Match_in_APC_functional_df == TRUE)

all_with_apc_data_clean <- all_with_apc_data[, c(2, 36:69)]

## Subset where there IS a Z-score (either from df3 OR df4).
all_with_z_score <- all_with_apc_data_clean %>%
  filter(!is.na(Z.score..df3.) | !is.na(Z.score..df4.))

## Subset where there IS NO Z-score (both df3 AND df4 are NA).
all_no_z_score <- all_with_apc_data_clean %>%
  filter(is.na(Z.score..df3.) & is.na(Z.score..df4.))

### Determine functional data evidence strengths (PS3/BS3) based on Thompson et al. (2024).
###########################################################################################
## Define a dataframe to store the classification results for variants with a Z score.
Functional_evidence_result <- data.frame(Final_cnomen = character(), 
                                         Final_gnomen = character(), 
                                         Final_pnomen = character(),
                                         Z_score_df4 = numeric(),
                                         Z_score_df3 = numeric(),
                                         Z_score_used = numeric(),
                                         PS3_Strong = logical(), 
                                         PS3_Mod = logical(), 
                                         PS3_Supp = logical(), 
                                         BS3_Supp = logical(), 
                                         BS3_Mod = logical())

## Loop through all variants with a Z score and determine criteria strengths based on thresholds.
for (i in seq_len(nrow(all_with_z_score))) {
  
  Variant_row <- all_with_z_score[i, ]
  
  Final_cnomen <- Variant_row$Final_cnomen
  Final_gnomen <- Variant_row$Final_gnomen
  Final_pnomen <- Variant_row$Final_pnomen
  Z_score_df4  <- Variant_row$Z.score..df4.
  Z_score_df3  <- Variant_row$Z.score..df3.
  
  # Select Z-score from df4 preferentially; if not available, use df3
  if (!is.na(Z_score_df4)) {
    Z_score_used <- Z_score_df4
  } else {
    Z_score_used <- Z_score_df3
  }
  
  # Initialize logical flags for criteria levels
  PS3_Strong <- FALSE
  PS3_Mod    <- FALSE
  PS3_Supp   <- FALSE
  BS3_Supp   <- FALSE
  BS3_Mod    <- FALSE
  
  # Categorize evidence strength according to the robust patch-clamp scale
  if (Z_score_used < -4) {
    PS3_Strong <- TRUE
  } else if (Z_score_used >= -4 & Z_score_used < -3) {
    PS3_Mod    <- TRUE
  } else if (Z_score_used >= -3 & Z_score_used <= -2) {
    PS3_Supp   <- TRUE
  } else if (Z_score_used > 1 | (Z_score_used >= -2 & Z_score_used < -1)) {
    BS3_Supp   <- TRUE
  } else if (Z_score_used >= -1 & Z_score_used <= 1) {
    BS3_Mod    <- TRUE
  }
  
  # Create a final result row and append to the cumulative dataframe
  final_row <- data.frame(Final_cnomen = Final_cnomen, 
                          Final_gnomen = Final_gnomen, 
                          Final_pnomen = Final_pnomen,
                          Z_score_df4 = Z_score_df4,
                          Z_score_df3 = Z_score_df3,
                          Z_score_used = Z_score_used,
                          PS3_Strong = PS3_Strong, 
                          PS3_Mod = PS3_Mod, 
                          PS3_Supp = PS3_Supp, 
                          BS3_Supp = BS3_Supp, 
                          BS3_Mod = BS3_Mod)
  
  Functional_evidence_result <- rbind(Functional_evidence_result, final_row)
}

## Save this result.
write.csv(Functional_evidence_result, "Functional_evidence_result_KCNH2_PS3_BS3.csv", row.names = FALSE, quote = TRUE)
