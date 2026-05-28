### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load in required libraries.
###############################
library(dplyr)
library(stringr)

### Load in variant data. 
#############################
Variants_combined_with_APC_functional_data_unique <- read.csv("Variants_combined_with_APC_functional_data_unique.csv")
Results_variant_recoder_APC_Z_filtered <- read.csv("Results_variant_recoder_APC_Z_filtered.csv")

## Ensure variant list column names are standardized for string matching
Variants_combined_with_APC_functional_data_unique$Final_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", Variants_combined_with_APC_functional_data_unique$Final_cnomen)

### Identify and filter APC variants. 
####################################
## Keep all variants where a match exists in the APC functional dataframe (irrespective of classification)
all_unique_base <- Variants_combined_with_APC_functional_data_unique %>%
  filter(Match_in_APC_functional_df == TRUE)

## Keep only variants where there is a valid functional Z-score
all_unique_base <- all_unique_base %>%
  filter(!is.na(Z.score..df3.) | !is.na(Z.score..df4.))


### Load the edited Oxford classification document.
###################################################
Oxford_classifications_edited <- read.csv("Oxford_classifications_for_all_lab_variants - EDIT.csv")

## Filter out entries with invalid placeholder labels in the Final classification column
Oxford_classifications_edited <- Oxford_classifications_edited %>%
  filter(
    !is.na(Final_classification) &
      trimws(Final_classification) != "" &
      trimws(Final_classification) != "N/A" &
      trimws(Final_classification) != "N/A."
  )

## Filter the Oxford document to keep all corresponding functional variants
Oxford_edited_filtered <- Oxford_classifications_edited %>%
  filter(Variant_cnomen %in% all_unique_base$Final_cnomen)


### Summarise classifications from the updated variant recoder dataframe.
##########################################################################
## Standardize cDNA column nomenclature by stripping prefixes to prevent mismatching
Results_variant_recoder_APC_Z_filtered$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", Results_variant_recoder_APC_Z_filtered$Final_cnomen)

## Extract and collapse unique classifications per variant across all lab submissions
Lab_classifications_collapsed <- Results_variant_recoder_APC_Z_filtered %>%
  group_by(Match_cnomen) %>%
  summarise(
    Collated_Classification_String = {
      vals <- trimws(na.omit(Lab.variant.classification))
      vals <- vals[vals != ""]
      if (length(vals) == 0) "" else paste(unique(vals), collapse = " | ")
    },
    .groups = "drop"
  )


### Load in all individual criteria data frames.
#################################################
Functional_evidence_result <- read.csv("Functional_evidence_result_KCNH2_PS3_BS3.csv")
Allele_freq_criterion_result <- read.csv("Allele_freq_criterion_result_GnomAD_260426.csv")
REVEL_evidence_result <- read.csv("REVEL_evidence_result_KCNH2.csv")
SpliceAI_evidence_result <- read.csv("SpliceAI_evidence_result_KCNH2.csv")
PM1_domain_criterion_results <- read.csv("KCNH2_PM1_hotspot_evaluations.csv")

## Clean up merge keys by stripping accession numbers across datasets
Functional_evidence_result$Match_key   <- gsub("^NM_000238(\\.[0-9]+)?:", "", Functional_evidence_result$Final_cnomen)
Allele_freq_criterion_result$Match_key <- gsub("^NM_000238(\\.[0-9]+)?:", "", Allele_freq_criterion_result$HGVSC)
REVEL_evidence_result$Match_key        <- gsub("^NM_000238(\\.[0-9]+)?:", "", REVEL_evidence_result$Variant_cnomen)
SpliceAI_evidence_result$Match_key     <- gsub("^NM_000238(\\.[0-9]+)?:", "", SpliceAI_evidence_result$Variant_cnomen)
PM1_domain_criterion_results$Match_key <- gsub("^NM_000238(\\.[0-9]+)?:", "", PM1_domain_criterion_results$Final_cnomen)


### Initialize the final compilation matrix dataframe layout.
###############################################################
## Structural Specifications layout rearranged to group variant tracking keys at the start
Combined_criteria_matrix_results <- data.frame(Variant_cnomen = character(), 
                                               Final_pnomen = character(),
                                               Collated_classification = character(), # Moved to follow final_pnomen
                                               PVS1 = character(), PS1 = character(), 
                                               PS2 = character(), PS3 = character(), PS4 = character(), 
                                               PM1 = character(), PM4 = character(), PM5 = character(), 
                                               PM6 = character(), PP1 = character(), PP2 = character(), 
                                               PP3 = character(), PM2 = character(), BA1 = character(), 
                                               BS1 = character(), BS3 = character(), BS4 = character(), 
                                               BP2 = character(), BP4 = character(), BP7 = character(),
                                               Total_ACGS_Points = numeric(),
                                               Final_classification = character(),    
                                               Comment = character(),
                                               stringsAsFactors = FALSE)


### Loop through all identified functional variants with Z-scores to compile criteria values.
#############################################################################################
for (i in seq_len(nrow(all_unique_base))) {
  
  Base_row  <- all_unique_base[i, ]
  c_nomen   <- Base_row$Final_cnomen
  raw_cname <- Base_row$Final_cnomen
  raw_pname <- Base_row$Final_pnomen
  
  # Pull criteria rows matching current variant key
  Func_row <- Functional_evidence_result[Functional_evidence_result$Match_key == c_nomen, ]
  Freq_row <- Allele_freq_criterion_result[Allele_freq_criterion_result$Match_key == c_nomen, ]
  Rev_row  <- REVEL_evidence_result[REVEL_evidence_result$Match_key == c_nomen, ]
  Spl_row  <- SpliceAI_evidence_result[SpliceAI_evidence_result$Match_key == c_nomen, ]
  PM1_row  <- PM1_domain_criterion_results[PM1_domain_criterion_results$Match_key == c_nomen, ]
  
  # Initialize empty character string cells for output matrix matching the template
  PVS1_val <- ""; PS1_val <- ""; PS2_val <- ""; PS3_val <- ""; PS4_val <- ""
  PM1_val  <- ""; PM4_val <- ""; PM5_val <- ""; PM6_val <- ""; PP1_val  <- ""
  PP2_val  <- ""; PP3_val  <- ""; PM2_val  <- ""; BA1_val  <- ""; BS1_val  <- ""
  BS3_val  <- ""; BS4_val  <- ""; BP2_val  <- ""; BP4_val  <- ""; BP7_val  <- ""
  
  # Initialize empty string tracking variables
  Final_class_val <- ""
  Collated_class_val <- ""
  
  # Pull the collated classifications associated with this variant if available
  Match_class_row <- Lab_classifications_collapsed[Lab_classifications_collapsed$Match_cnomen == c_nomen, ]
  if (nrow(Match_class_row) > 0) {
    Collated_class_val <- Match_class_row$Collated_Classification_String[1]
  }
  
  # Initialize point counter variable
  Points <- 0
  
  # 1. Translate PS3 and BS3 strengths from Functional data loop outputs
  if (nrow(Func_row) > 0) {
    if (any(Func_row$PS3_Strong, na.rm = TRUE)) { PS3_val <- "Str"; Points <- Points + 4 }
    if (any(Func_row$PS3_Mod, na.rm = TRUE))    { PS3_val <- "Mod"; Points <- Points + 2 }
    if (any(Func_row$PS3_Supp, na.rm = TRUE))   { PS3_val <- "Sup"; Points <- Points + 1 }
    if (any(Func_row$BS3_Supp, na.rm = TRUE))   { BS3_val <- "Sup"; Points <- Points - 1 }
    if (any(Func_row$BS3_Mod, na.rm = TRUE))    { BS3_val <- "Mod"; Points <- Points - 2 }
  }
  
  # 2. Translate PM2, BS1, BA1 from Allele Frequency loop outputs
  if (nrow(Freq_row) > 0) {
    if (any(Freq_row$PM2_Supporting, na.rm = TRUE)) { PM2_val <- "Sup"; Points <- Points + 1 }
    if (any(Freq_row$BS1, na.rm = TRUE))            { BS1_val <- "Str"; Points <- Points - 4 }
    if (any(Freq_row$BA1, na.rm = TRUE))            { BA1_val <- "VStr"; Points <- Points - 8 } 
  }
  
  # 3. Translate PM1 Hotspot strings from Domain loop outputs
  if (nrow(PM1_row) > 0) {
    if (any(PM1_row$PM1_mod, na.rm = TRUE)) { PM1_val <- "Mod"; Points <- Points + 2 }
    if (any(PM1_row$PM1_sup, na.rm = TRUE)) { PM1_val <- "Sup"; Points <- Points + 1 }
  }
  
  # 4. Translate REVEL and SpliceAI into computational cells PP3 / BP4
  Comp_sources <- c()
  
  if (nrow(Rev_row) > 0) {
    if (any(Rev_row$PP3, na.rm = TRUE)) { 
      PP3_val <- "Sup" 
      Comp_sources <- c(Comp_sources, "REVEL (PP3)")
    }
    if (any(Rev_row$BP4, na.rm = TRUE)) { 
      BP4_val <- "Sup" 
      Comp_sources <- c(Comp_sources, "REVEL (BP4)")
    }
  }
  
  if (nrow(Spl_row) > 0) {
    if (any(Spl_row$PP3, na.rm = TRUE)) { 
      PP3_val <- "Sup" 
      Comp_sources <- c(Comp_sources, "SpliceAI (PP3)")
    }
  }
  
  # Add points once if either/both computational prediction platforms trigger a rule
  if (PP3_val == "Sup") { Points <- Points + 1 }
  if (BP4_val == "Sup") { Points <- Points - 1 }
  
  # Initialize comment value or format the vector into a text comment
  Comment_val <- ""
  if (length(Comp_sources) > 0) {
    Comment_val <- paste0("Computational data applied: ", paste(Comp_sources, collapse = ", "))
  }
  
  # Compile row matching template columns matching structural specifications ordering
  final_matrix_row <- data.frame(Variant_cnomen = raw_cname, 
                                 Final_pnomen = raw_pname,
                                 Collated_classification = Collated_class_val, # Position adjusted
                                 PVS1 = PVS1_val, PS1 = PS1_val, 
                                 PS2 = PS2_val, PS3 = PS3_val, PS4 = PS4_val, 
                                 PM1 = PM1_val, PM4 = PM4_val, PM5 = PM5_val, 
                                 PM6 = PM6_val, PP1 = PP1_val, PP2 = PP2_val, 
                                 PP3 = PP3_val, PM2 = PM2_val, BA1 = BA1_val, 
                                 BS1 = BS1_val, BS3 = BS3_val, BS4 = BS4_val, 
                                 BP2 = BP2_val, BP4 = BP4_val, BP7 = BP7_val,
                                 Total_ACGS_Points = Points,
                                 Final_classification = Final_class_val,       
                                 Comment = Comment_val,
                                 stringsAsFactors = FALSE)
  
  Combined_criteria_matrix_results <- rbind(Combined_criteria_matrix_results, final_matrix_row)
}

## Save this result.
write.csv(Combined_criteria_matrix_results, "Combined_criteria_matrix_results.csv", row.names = FALSE, quote = TRUE)


### Filter the final variant recoder data frame by the identified cohort variants.
#################################################################################

## Filter to only keep rows matching the unique variants with valid functional data from the cohort
Results_variant_recoder_APC_Z_filtered <- Results_variant_recoder_APC_Z_filtered %>%
  filter(Match_cnomen %in% all_unique_base$Final_cnomen)

## Clean up by removing the temporary matching column
Results_variant_recoder_APC_Z_filtered$Match_cnomen <- NULL

## Save this final subset dataframe result.
write.csv(Results_variant_recoder_APC_Z_filtered, "Results_variant_recoder_APC_Z_filtered.csv", row.names = FALSE, quote = TRUE)


### Filter the Oxford classifications dataframe by useable ACGS criteria records.
##############################################################################

## Filter the filtered cohort dataset for rows containing useable classification criteria
useable_criteria_df <- Results_variant_recoder_APC_Z_filtered %>%
  filter(
    !is.na(ACGS.criteria) & 
      trimws(ACGS.criteria) != "" & 
      trimws(ACGS.criteria) != "N/A" & 
      trimws(ACGS.criteria) != "Pre ACMG guidelines"
  )

## Standardize cDNA column nomenclature by stripping prefixes to prevent mismatching
useable_criteria_df$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", useable_criteria_df$Final_cnomen)

## Filter the raw Oxford classifications document to keep only these matching records
Oxford_filtered_by_ACGS_criteria <- Oxford_classifications_edited %>%
  filter(Variant_cnomen %in% useable_criteria_df$Match_cnomen)

## Save this specific filtered result.
write.csv(Oxford_filtered_by_ACGS_criteria, "Oxford_filtered_by_ACGS_criteria.csv", row.names = FALSE, quote = TRUE)

## Save this result.
# NOTE. Variant c.388G>A is written as VUS when upon further inspection this was recorded as likely pathogenic. Edit this in VUS list. 
write.csv(Combined_criteria_matrix_results, "Combined_criteria_matrix_results.csv", row.names = FALSE, quote = TRUE)
## Edit this to assess the VUS variants and their classification in Oxford and external labs
## to determine if any extra criteria could be applied. 

### Filter the final variant recoder data frame by the identified cohort variants.
#################################################################################

## Filter to only keep rows matching the unique variants with valid functional data from the cohort
Results_variant_recoder_APC_Z_filtered <- Results_variant_recoder_protein_df_updated %>%
  filter(Match_cnomen %in% all_unique_base$Final_cnomen)

## Clean up by removing the temporary matching column
Results_variant_recoder_APC_Z_filtered$Match_cnomen <- NULL

## Save this final subset dataframe result.
write.csv(Results_variant_recoder_APC_Z_filtered, "Results_variant_recoder_APC_Z_filtered.csv", row.names = FALSE, quote = TRUE)


### Filter the Oxford classifications dataframe by useable ACGS criteria records.
##############################################################################

## Filter the filtered cohort dataset for rows containing useable classification criteria
useable_criteria_df <- Results_variant_recoder_APC_Z_filtered %>%
  filter(
    !is.na(ACGS.criteria) & 
      trimws(ACGS.criteria) != "" & 
      trimws(ACGS.criteria) != "N/A" & 
      trimws(ACGS.criteria) != "Pre ACMG guidelines"
  )

## Standardize cDNA column nomenclature by stripping prefixes to prevent mismatching
useable_criteria_df$Match_cnomen <- gsub("^NM_000238(\\.[0-9]+)?:", "", useable_criteria_df$Final_cnomen)

## Filter the raw Oxford classifications document to keep only these matching records
Oxford_filtered_by_ACGS_criteria <- Oxford_classifications_edited %>%
  filter(Variant_cnomen %in% useable_criteria_df$Match_cnomen)

## Save this filtered result.
write.csv(Oxford_filtered_by_ACGS_criteria, "Oxford_filtered_by_ACGS_criteria.csv", row.names = FALSE, quote = TRUE)