### Script to create collated variant dataset from data recieved from Aberdeen, Brompton, 
### Belfast, Manchester, and own data at Oxford. 

### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")
Projectdir <- "C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/"

### Load in required libraries.
###############################
library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble) 

### Load in lab variant data.
#############################

## Aberdeen. Using updated spreadsheet sent via email along with variant correction described below. 
Aberdeen <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/Aberdeen - KCNH2 variants - Updated_May 2025.xlsx", sep = ""), skip = 2)
Aberdeen <- Aberdeen[-c(1,29),]

# Assign lab name. 
Aberdeen$Lab <- "Aberdeen"

# Received correspondence from Aberdeen re: variant c.65C>T. This was corrected to c.65T>C. Change this. 
Aberdeen[43,]$`Variant cDNA no` <- "c.65T>C"

# Fix column names to easily merge with other lab variant lists. 
colnames(Aberdeen) <- c("c.nomenclature", "p.nomenclature", "Lab variant classification", "Classification date", 
                        "Occurrences", "Segregation/ Meiosis  (excluding proband)", "Panel", "ACGS criteria",
                        "Lab")

## Brompton. 
# VUS variants. 
Brompton_VUS <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/Brompton - KCNH2 variants - RBH.xlsx", sep = ""), sheet = "Class - VUS")
Brompton_VUS$Classification <- "VUS"
colnames(Brompton_VUS) <- c("LabNo", "notes", "c.", "p.", "finalReportAuthorized", "panelReported", "Classification")

# Likely pathogenic variants.
Brompton_LikelyPath <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/Brompton - KCNH2 variants - RBH.xlsx", sep = ""), sheet = "Class - Likely Pathogenic")
Brompton_LikelyPath$Classification <- "Likely pathogenic"
colnames(Brompton_LikelyPath) <- c("LabNo", "notes", "c.", "p.", "finalReportAuthorized", "panelReported", "Classification")

# Pathogenic variants. 
Brompton_Path <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/Brompton - KCNH2 variants - RBH.xlsx", sep = ""), sheet = "Class -  Pathogenic")
Brompton_Path <- Brompton_Path[,-c(7:10)]
Brompton_Path$Classification <- "Pathogenic"
colnames(Brompton_Path) <- c("LabNo", "notes", "c.", "p.", "finalReportAuthorized", "panelReported", "Classification")

# Bind rows together for each classification. 
Brompton <- bind_rows(Brompton_VUS, Brompton_LikelyPath, Brompton_Path)

# Change date from US to UK format. 
Brompton[[5]] <- format(Brompton[[5]], "%d-%m-%Y")

# Remove leading and trailing commas. 
Brompton$c. <- gsub("^[.,]+|[.,]+$", "", Brompton$c.)
Brompton$p. <- gsub("^[.,]+|[.,]+$", "", Brompton$p.)

# Split entries where two variants were detected. 
Brompton <- Brompton %>%
  mutate(
    c_list = str_split(`c.`, "[,;]"),
    p_list = str_split(`p.`, "[,;]")
  ) %>%
  mutate(var_pairs = map2(c_list, p_list, ~ tibble(`c.` = .x, `p.` = .y))) %>%
  select(LabNo, notes, everything(), -`c.`, -`p.`, -c_list, -p_list) %>%
  unnest(var_pairs) %>%
  group_by(LabNo) %>%
  mutate(
    other_c = map_chr(`c.`, ~ paste(setdiff(`c.`, .x), collapse = ",")),
    notes = case_when(
      other_c != "" ~ if_else(
        is.na(notes) | notes == "",
        paste0("Detected alongside ", other_c),
        paste0(notes, " | Detected alongside ", other_c)
      ),
      TRUE ~ notes
    )
  ) %>%
  ungroup() %>%
  select(-other_c)

# Received correspondence from Brompton that c.233A>G variant is in KCNJ2, so not a ?comp het KCNH2 case. 
# Delete the KCNH2 row for c.233A>G and alter note for c.2362G>A variant. 
Brompton <- Brompton[-80,]
Brompton[79,]$notes <- "Detected alongside 	KCNJ2 c.233A>G p.(Asp78Gly)"

# Edits for further correspondence which provided more information on the minimally deleted regions according 
# to the used CNV tool in Brompton for deletions in variant list provided. 
Brompton[133,]$notes <- "Minimally deleted region chr7:150946876-150947514. Proband, child of LabNo"
Brompton[150,]$notes <- "Minimally deleted region chr7:150946876-150947514. Parent of LabNo"

Brompton[155,]$notes <- "Minimally deleted region chr7:150945364-150945514."

# Other edits for corrected nomenclature.
Brompton[30,]$p. <- "p.(Arg176Trp)"

Brompton[4,]$p. <- "p.(Ala913Val)"

Brompton$p.[Brompton$p. == "p.Asp629Ser"] <- "p.Asn629Ser"

Brompton$c.[Brompton$c. == "c.3090_3012del"] <- "c.3090_3102del"

# Assign lab name. 
Brompton$Lab <- "Brompton"

# Fix column names to easily merge with other lab variant lists. 
colnames(Brompton) <- c("Lab number", "Notes", "Classification date", "Panel", "Lab variant classification", 
                        "c.nomenclature", "p.nomenclature", "Lab")

## Belfast. 
Belfast <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/KCNH2_Belfast - Updated_Aug 2025.xlsx", sep = ""), col_types = c("text", "text", "text", "text", "text", 
                                                                                                                                               "text","date", "text", "text", 
                                                                                                                                               "text", "text", "text"))

# Change date formatting from US to UK. 
Belfast[[7]] <- format(Belfast[[7]], "%d-%m-%Y")
Belfast[c(12, 14), 7] <- as.character("2008")

# Remove empty row. 
Belfast <- Belfast[-23,]

# Assign lab name. 
Belfast$Lab <- "Belfast"

# Fix column names to easily merge with other lab variant lists. 
colnames(Belfast) <- c("Gene", "Transcript", "c.nomenclature", "p.nomenclature", "Lab variant classification", 
                        "ACGS criteria", "Classification date", "Occurrences", "Panel", "Family members with variant",
                       "Family members with known phenotype", "Notes", "Lab")

## Manchester. 
Manchester <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/Manchester_KCNH2 Data Request DRQ 00091.xlsx", sep = ""), col_names = FALSE)

# Variants sent in single column as strings. Split these into appropriate columns. Variants detected 
# alongside the KCNH2 variant have also been noted. These are noted in the comments column. 
Manchester <- Manchester %>%
  mutate(all_variants = str_split(...1, ",\\s*")) %>%       # Split variants by comma
  mutate(
    KCNH2_variants = map(all_variants, ~ .x[str_detect(.x, "^KCNH2")]),  # Extract KCNH2
    other_variants = map(all_variants, ~ .x[!str_detect(.x, "^KCNH2")])  # Extract others
  ) %>%
  unnest_longer(KCNH2_variants, values_to = "KCNH2_var", keep_empty = TRUE) %>% # One row per KCNH2 variant (or NA if none)
  mutate(
    Gene = if_else(!is.na(KCNH2_var), "KCNH2", NA_character_),
    c. = str_extract(KCNH2_var, "c\\.[^ ]+"),
    p. = str_extract(KCNH2_var, "p\\.\\([^\\)]+\\)"),
    Zygosity = str_extract(KCNH2_var, "(?i)het|hom"),
    comments = map_chr(other_variants, ~ paste(.x, collapse = ", "))   # Combine non-KCNH2 variants into comments
  ) %>%
  select(-all_variants, -other_variants, -KCNH2_var)

Manchester <- Manchester %>%
  mutate(
    comments = if_else(
      comments != "" & !is.na(comments),
      paste0("Detected alongside ", comments),
      comments
    )
  )

# For entries where there are two KCNH2 variants detected, a row has been created for each unique variant, but the 
# additional KCNH2 variant has been noted in the comments tab. In each instance this occurs, the patient has been 
# assigned a number starting from one so variants detected together in an individual stay linked. 
Manchester[33,]$comments <- c("Detected alongside KCNH2 c.2230C>T p.(Arg744Ter) het? Assigning as patient 1")
Manchester[34,]$comments <- c("Detected alongside KCNH2 c.2738C>T p.(Ala913Val) het? Assigning as patient 1")
Manchester[61,]$comments <- c("Detected alongside KCNH2 c.2690A>C p.(Lys897Thr) het? Assigning as patient 2")
Manchester[62,]$comments <- c("Detected alongside KCNH2 c.2860C>T p.(Arg954Cys) het? Assigning as patient 2")
Manchester[64,]$comments <- c("Detected alongside KCNH2 c.2582A>T p.(Asn861Ile) het? Assigning as patient 3")
Manchester[65,]$comments <- c("Detected alongside KCNH2 c.961G>C p.(Asp321His) het Assigning as patient 3")

# Assign lab name. 
Manchester$Lab <- "Manchester"

# Fix column names to easily merge with other lab variant lists. 
colnames(Manchester) <- c("Original entry", "Gene", "c.nomenclature", "p.nomenclature", "Zygosity", "Notes",
                          "Lab")

## Oxford.
Oxford <- read_xlsx(paste(Projectdir,"Variants/Variants from labs (unedited)/KCNH2_ORGL_2021.xlsx", sep = ""))

# Assign lab name. 
Oxford$Lab <- "Oxford"

# Fix column names to easily merge with other lab variant lists. 
colnames(Oxford) <- c("Gene", "g.nomenclature", "c.nomenclature", "p.nomenclature", "Lab variant classification", 
                      "Occurrences", "Unique?", "Updated pathogenicity", "KT checked", "KT-class with new PM1 data",
                      "KT Comment", "Lab")

# Replace classification column with 'updated classification'. This will still need re-retrieving later. 
Oxford$`Lab variant classification`[!is.na(Oxford$`Updated pathogenicity`) & Oxford$`Updated pathogenicity` != ""] <- 
  Oxford$`Updated pathogenicity`[!is.na(Oxford$`Updated pathogenicity`) & Oxford$`Updated pathogenicity` != ""]

# Edit for corrected nomenclature.
Oxford$c.nomenclature[Oxford$c.nomenclature == "c.982_983ins13"] <- "c.982_983insATTATGCGCTAC"

### Merge all variants from different labs. 
###########################################

## Create empty data frames containing all relevant columns needed across the different 
## lab variant lists. 
Collated_variant_list <- data.frame(
  `Lab` = character(),
  `Lab number` = character(),
  `Panel` = character(),
  `Transcript` = character(),
  `g.nomenclature` = character(),
  `c.nomenclature` = character(),
  `p.nomenclature` = character(),
  `Lab variant classification` = character(),
  `Occurrences` = character(),
  `Unique?` = character(),
  `Classification date` = character(),
  `Oxford current classification` = character(),
  `Functional data available?` = character(),
  `ACGS criteria` = character(),
  `Family members with variant` = character(),
  `Family members with known phenotype` = character(),
  `Notes` = character(),
  stringsAsFactors = FALSE, check.names = FALSE
)

## Merge the data frames. 
# Create a list of the data frames to merge. 
dfs <- list(Collated_variant_list, Aberdeen, Belfast, Brompton, Manchester, Oxford)

# Create character list of all the columns to include in final data frame. 
template_cols <- names(Collated_variant_list)

# Make each data frame match the template (Collated_variant_list) (add missing cols as 
# NA, drop extras, order columns)
dfs_trimmed <- lapply(dfs, function(x) {
  missing <- setdiff(template_cols, names(x))
  if (length(missing)) x[missing] <- NA_character_
  x <- x[template_cols]   
  return(x)
})

# Merge them into one big df
Collated_variant_list <- do.call(rbind, dfs_trimmed)

## Normalise naming of classifications.
Collated_variant_list <- Collated_variant_list %>%
  mutate(`Lab variant classification` = case_when(
    `Lab variant classification` %in% c("Class 3 (VUS)", "VUS", "Class 3-Unknown pathogenicity", "Class 3-Unknown pathogenicity (hot)", "Class 3-Unknown pathogenicity (cold)") ~ "VUS",
    `Lab variant classification` %in% c("Class 4", "Likely pathogenic", "LP", "Class 4-Likely pathogenic") ~ "Likely pathogenic",
    `Lab variant classification` %in% c("P", "Pathogenic", "Class 5-Certainly pathogenic") ~ "Pathogenic",
    `Lab variant classification` %in% c("LB") ~ "Likely benign",
    `Lab variant classification` %in% c("B") ~ "Benign",
    `Lab variant classification` %in% c("Risk alele") ~ "Risk allele",
    TRUE ~ as.character(`Lab variant classification`)  # keep other values as-is, including NA
  ))

## Remove brackets from p.nomenclature.
Collated_variant_list$p.nomenclature <- gsub("[()]", "", Collated_variant_list$p.nomenclature)

## Fix entries misspelled. 
Collated_variant_list$p.nomenclature[Collated_variant_list$p.nomenclature == "p.lle789Phe"] <- "p.Ile789Phe"
Collated_variant_list$p.nomenclature[Collated_variant_list$p.nomenclature == "p.lle19Thr"] <- "p.Ile19Thr"
Collated_variant_list$p.nomenclature[Collated_variant_list$p.nomenclature == "p.lle82_Gln84dup"] <- "p.Ile82_Gln84dup"

Collated_variant_list[is.na(Collated_variant_list)] <- ""

## Save variant list. 
write.csv(Collated_variant_list, "Collated_variant_list.csv", row.names = FALSE, quote = TRUE)
