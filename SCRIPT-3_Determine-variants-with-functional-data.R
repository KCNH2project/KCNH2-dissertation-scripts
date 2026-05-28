### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")
Projectdir <- "C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/"
Functional_data_dir <- "C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/Data from functional studies"

### Load in required libraries.
###############################
library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble) 


### Load in functional data.
############################

## A calibrated functional patch-clamp assay to enhance clinical variant interpretation in KCNH2-related long QT syndrome, Connie Jiang et al (2022). Patch clamp. 
A_calibrated_functional_Connie_Jiang_et_al_2022_VUS <- read_xlsx(paste(Functional_data_dir,"/A_calibrated_functional-Connie_Jiang_et_al_2022.xlsx", 
                                                                     sep = ""), sheet = "Table S4_EP data", range="B42:G88")

A_calibrated_functional_Connie_Jiang_et_al_2022_Benign <- read_xlsx(paste(Functional_data_dir,"/A_calibrated_functional-Connie_Jiang_et_al_2022.xlsx", 
                                                                     sep = ""), sheet = "Table S4_EP data", range="B3:G18")

A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic <- read_xlsx(paste(Functional_data_dir,"/A_calibrated_functional-Connie_Jiang_et_al_2022.xlsx", 
                                                                     sep = ""), sheet = "Table S4_EP data", range="B21:G39")

## A massively parallel assay accurately discriminates between functionally normal and abnormal variants in a hotspot domain of KCNH2, Chai-Ann Ng et al (2022). Patch clamp.
A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022 <- read_xlsx(paste(Functional_data_dir,"/A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022.xlsx", 
                                                                      sep = ""), sheet = "Table S1", range="B2:L462", .name_repair = "minimal")

## Clinical interpretation of KCNH2 variants using a robust PS3/BS3 functional patch-clamp assay, Kate L Thomson et al (2024). Patch clamp.
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS <- read_xlsx(paste(Functional_data_dir,"/Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024.xlsx", 
                                                                                sep = ""), sheet = "Table S3_Patch_Clamp_data", range="B73:G92", .name_repair = "minimal")
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign <- read_xlsx(paste(Functional_data_dir,"/Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024.xlsx", 
                                                                                sep = ""), sheet = "Table S3_Patch_Clamp_data", range="B3:G34", .name_repair = "minimal")
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic <- read_xlsx(paste(Functional_data_dir,"/Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024.xlsx", 
                                                                                sep = ""), sheet = "Table S3_Patch_Clamp_data", range="B38:G69", .name_repair = "minimal")

## Multiplexed Assays of Variant Effect and Automated Patch Clamping Improve KCNH2-LQTS Variant Classification and Cardiac Event Risk Stratification, Matthew J O'Neill et al (2024). Patch clamp. 
Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC <- read_xlsx(paste(Functional_data_dir,"/Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024.xlsx", 
                                                                              sep = ""), sheet = "Table S5", range="A3:I537", .name_repair = "minimal")

### Tidy functional datasets.
#############################

## Tidying A_calibrated... 
# Clean A_calibrated_functional_Connie_Jiang_et_al_2022_VUS.
colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_VUS) <- paste0(colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_VUS), 
                                                                      " ", as.character(A_calibrated_functional_Connie_Jiang_et_al_2022_VUS[1, ]))
colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_VUS) <- gsub(" NA", "", colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_VUS))
A_calibrated_functional_Connie_Jiang_et_al_2022_VUS <- A_calibrated_functional_Connie_Jiang_et_al_2022_VUS[-1,]
A_calibrated_functional_Connie_Jiang_et_al_2022_VUS$`Protein Change` <- gsub("[*^]", "", A_calibrated_functional_Connie_Jiang_et_al_2022_VUS$`Protein Change`)

# Add column which highlights this is a list of VUS in the study. 
A_calibrated_functional_Connie_Jiang_et_al_2022_VUS$Study_classification <- "VUS to be investigated."

# Same for A_calibrated_functional_Connie_Jiang_et_al_2022_Benign.
colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Benign) <- paste0(colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Benign), 
                                                                      " ", as.character(A_calibrated_functional_Connie_Jiang_et_al_2022_Benign[1, ]))
colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Benign) <- gsub(" NA", "", colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Benign))
A_calibrated_functional_Connie_Jiang_et_al_2022_Benign <- A_calibrated_functional_Connie_Jiang_et_al_2022_Benign[-1,]

# Add column which highlights this is a list of benign controls in the study. 
A_calibrated_functional_Connie_Jiang_et_al_2022_Benign$Study_classification <- "Benign control."

# Same for A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic
colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic) <- paste0(colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic), 
                                                                      " ", as.character(A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic[1, ]))
colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic) <- gsub(" NA", "", colnames(A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic))
A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic <- A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic[-1,]

# Add column which highlights this is a list of pathogenic controls in the study. 
A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic$Study_classification <- "Pathogenic control."

# Merge A_calibrated_functional_Connie_Jiang_et_al_2022_VUS, A_calibrated_functional_Connie_Jiang_et_al_2022_Benign, and A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic
# protein change. 
dfs <- list(A_calibrated_functional_Connie_Jiang_et_al_2022_VUS, A_calibrated_functional_Connie_Jiang_et_al_2022_Benign, 
            A_calibrated_functional_Connie_Jiang_et_al_2022_Pathogenic)

# Merge on all shared columns (e.g., "Protein Change", etc.)
shared_cols <- Reduce(intersect, lapply(dfs, names))

A_calibrated_functional_Connie_Jiang_et_al_2022 <- reduce(dfs, full_join, by = shared_cols)

# Add column containing study name. 
A_calibrated_functional_Connie_Jiang_et_al_2022$'Study name' <- "Jiang, C. et al. (2022) “A calibrated functional patch-clamp assay to enhance clinical variant interpretation in KCNH2-related long QT syndrome,” American Journal of Human Genetics, 109(7), pp. 1199–1207."

## Tidying A_massively... 
# Fix column header names and remove the unnecessary rows.  
colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022) <- paste(colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022), 
                                                                      as.character(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022[1, ]), 
                                                                      as.character(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022[2, ]), 
                                                                      sep = " ")
colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022) <- gsub("\\s*\\bNA\\b\\s*", "", colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022))
A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022 <- A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022[-c(1,2),]
colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022)[5] <- "Variant Browser (Kozek et al., 2021) Affected LQTS"
colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022)[6] <- "Variant Browser (Kozek et al., 2021) Unaffected"
colnames(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022)[1] <- "Protein Change"

# Remove the p. at the beginning of the protein, to be added later for consistency. 
A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022$`Protein Change` <- gsub("[pP]\\.", "", A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022$`Protein Change`)

# Add column containing study name. 
A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022$'Study name' <- "Ng, C.A. et al. (2022) “A massively parallel assay accurately discriminates between functionally normal and abnormal variants in a hotspot domain of KCNH2,” American Journal of Human Genetics, 109(7), pp. 1208–1216."

## Tidying Clinical_interpretation...
# Fix column header names and remove the unnecessary row/column from Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign
colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign) <- paste0(colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign), 
                                                                                 " ", as.character(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign[1, ]))
colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign) <- gsub("\\s*\\bNA\\b\\s*", "", colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign))
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign <- Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign[-1,]

colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign)[6] <- "Functional Data Evidence Strength"

# Add column which highlights this is a list of benign controls in the study. 
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign$Study_classification <- "Benign control."

# Fix column header names and remove the unnecessary row/column from Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS
colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS) <- paste0(colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS), 
                                                                                 " ", as.character(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS[1, ]))
colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS) <- gsub("\\s*\\bNA\\b\\s*", "", colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS))
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS <- Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS[-1,]

colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS)[6] <- "Functional Data Evidence Strength"

# Add column which highlights this is a list of VUS in the study. 
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS$Study_classification <- "VUS to be investigated."

# Fix column header names and remove the unnecessary row/column from Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic
colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic) <- paste0(colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic), 
                                                                                 " ", as.character(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic[1, ]))
colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic) <- gsub("\\s*\\bNA\\b\\s*", "", colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic))
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic <- Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic[-1,]

colnames(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic)[6] <- "Functional Data Evidence Strength"

# Add column which highlights this is a list of pathogenic controls in the study. 
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic$Study_classification <- "Pathogenic control."

# Merge Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_1 and Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_2, _3, & _4 
# based on protein change. 
# Put your data frames in a list
dfs2 <- list(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Benign, Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_VUS, 
             Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_Pathogenic)

# Merge on all shared columns (e.g., "Protein Change", etc.)
shared_cols2 <- Reduce(intersect, lapply(dfs2, names))

Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024 <- reduce(dfs2, full_join, by = shared_cols2)

# Add column containing study name. 
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024$'Study name' <- "Thomson, K.L. et al. (2024) “Clinical interpretation of KCNH2 variants using a robust PS3/BS3 functional patch-clamp assay,” Human Genetics and Genomics Advances, 5(2), p. 100270."

# Remove the p. at the beginning of the protein, to be added later for consistency. 
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024$`Protein Change` <- gsub("[pP]\\.", "", Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024$`Protein Change`)

# Some rows are duplicated in this study as VUS were reclassified to become pathogenic or benign controls. Merge these into a single
# entry, keeping extra info in same cell separated by a ",". 
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024 <- Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024 %>%
  group_by(`Protein Change`) %>%
  summarise(across(
    everything(),
    ~ paste(unique(.), collapse = ", "),
    .names = "{.col}"
  ), .groups = "drop")

## Tidying Multiplexed_Assays...APC.
# Fix column header names and remove the unnecessary row/column.
colnames(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC) <- paste(colnames(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC), 
                                                                              as.character(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC[1, ]), 
                                                                              sep = " ")
colnames(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC) <- gsub("\\s*\\bNA\\b\\s*", "", colnames(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC))
Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC <- Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC[-1,-7]
colnames(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC)[1] <- "Protein Change"
colnames(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC)[8] <- "n number.1"

# Remove the p. at the beginning of the protein, to be added later for consistency. 
Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC$`Protein Change` <- gsub("[pP]\\.", "", Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC$`Protein Change`)

# Add column containing study name. 
Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC$'Study name' <- "O’Neill, M.J. et al. (2024) “Multiplexed Assays of Variant Effect and Automated Patch-clamping Improve KCNH2-LQTS Variant Classification and Cardiac Event Risk Stratification,” medRxiv, p. 2024.02.01.24301443."

### Merge the functional datasets.
##################################

# Helper function to suffix column names except the join column.
rename_with_suffix <- function(df, suffix) {
  names(df) <- ifelse(names(df) == "Protein Change", "Protein Change",
                      paste0(names(df), " ", suffix))
  return(df)
}

# Add suffixes
A_calibrated_functional_Connie_Jiang_et_al_2022_df1 <- rename_with_suffix(A_calibrated_functional_Connie_Jiang_et_al_2022, "(df1)")
A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022_df2 <- rename_with_suffix(A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022, "(df2)")
Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_df3 <- rename_with_suffix(Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024, "(df3)")
Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC_df4 <- rename_with_suffix(Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC, "(df4)")

# Merge all on 'Protein Change'
Merged_functional_df <- reduce(
  list(A_calibrated_functional_Connie_Jiang_et_al_2022_df1, 
       A_massively_parrallel_assay_Chai_Ann_Ng_et_al_2022_df2,
       Clinical_interpretation_of_KCNH2_Kate_L_Thomson_et_al_2024_df3, 
       Multiplexed_Assays_of_Variant_Effect_Matthew_J_ONeill_2024_APC_df4),
  function(x, y) full_join(x, y, by = "Protein Change")
)

# Combine study name columns
Merged_functional_df$'Study name' <- apply(
  Merged_functional_df[, grep("Study name*", names(Merged_functional_df))],
  1,
  function(x) paste(na.omit(unique(x)), collapse = " & ")
)

# Remove original study name columns
Merged_functional_df <- Merged_functional_df[, !grepl("^Study name", names(Merged_functional_df)) | names(Merged_functional_df) == "Study name"]

# Add p. to beginning of protein nomenclature column. 
Merged_functional_df$`Protein Change` <- paste0("p.", Merged_functional_df$`Protein Change`)

### Create csv from the merged data frame.
#write.csv(Merged_functional_df, paste0(getwd(),"/Merged_functional_data.csv"), row.names = FALSE)