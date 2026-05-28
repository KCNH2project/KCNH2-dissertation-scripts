### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load in required libraries.
###############################
library(stringr)

### Load in the variant data.
#############################
Results_variant_recoder_df <- read.csv("Results_variant_recoder_protein_df_final.csv")

### Remove duplicate variants based on cDNA nomenclature (Final_cnomen).
#######################################################################
Results_variant_recoder_df <- Results_variant_recoder_df[!duplicated(Results_variant_recoder_df$Final_cnomen), ]

### Initialize PM1 Domain Hotspot Dataframe.
############################################
PM1_domain_criterion_results <- data.frame(
  Final_gnomen = character(),
  Final_cnomen = character(),
  Final_pnomen = character(),
  PM1_mod = logical(),
  PM1_sup = logical(),
  Domain_Comment = character(),
  stringsAsFactors = FALSE
)

### Loop through variants to isolate codons and assign ACGS PM1 parameters.
###########################################################################
for (i in 1:nrow(Results_variant_recoder_df)) {
  
  # Extract variant identifiers
  Final_gnomen <- Results_variant_recoder_df$Final_gnomen[i]
  Final_cnomen <- Results_variant_recoder_df$Final_cnomen[i]
  Final_pnomen <- Results_variant_recoder_df$Final_pnomen[i]
  
  # Pull digits following the 'p.' prefix to determine the amino acid position
  Codon_extracted <- str_extract(Final_pnomen, "(?<=p\\.[A-Za-z]{3})\\d+|(?<=p\\.)\\d+")
  Codon <- as.numeric(Codon_extracted)
  
  # Default logic initialization (All variants outside specified domains default to FALSE / generic comments)
  PM1_mod        <- FALSE
  PM1_sup        <- FALSE
  Domain_Comment <- "Outside Hotspot Domains"
  
  # Handle exceptions where nomenclature is structural or missing (e.g. p.?)
  if (is.na(Codon)) {
    Domain_Comment <- "Non-missense / Missing nomenclature"
  } else {
    
    # Branch 1: N-terminus cluster (1 to 130) -> PM1 Moderate
    if (Codon >= 1 & Codon <= 130) {
      PM1_mod        <- TRUE
      Domain_Comment <- "N-terminus cluster"
      
      # Branch 2: Transmembrane/Linker/Pore (404 to 659) -> PM1 Moderate
    } else if (Codon >= 404 & Codon <= 659) {
      PM1_mod        <- TRUE
      Domain_Comment <- "Transmembrane/Linker/Pore"
      
      # Branch 3: C-terminus cNBD (742 to 842) -> PM1 Supporting
    } else if (Codon >= 742 & Codon <= 842) {
      PM1_sup        <- TRUE
      Domain_Comment <- "C-terminus cNBD"
    }
  }
  
  # Compile row entry matching exact layout request
  final_row <- data.frame(
    Final_gnomen   = Final_gnomen,
    Final_cnomen   = Final_cnomen,
    Final_pnomen   = Final_pnomen,
    PM1_mod        = PM1_mod,
    PM1_sup        = PM1_sup,
    Domain_Comment = Domain_Comment,
    stringsAsFactors = FALSE
  )
  
  # Bind back into final loop compilation matrix
  PM1_domain_criterion_results <- rbind(PM1_domain_criterion_results, final_row)
}

### Save the structured hotspot data.
#####################################
write.csv(PM1_domain_criterion_results, "KCNH2_PM1_hotspot_evaluations.csv", row.names = FALSE)


