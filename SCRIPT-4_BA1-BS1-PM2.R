### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load in required libraries.
###############################
library(tidyr)
library(dplyr)
library(stringr)
library(stats)
library(binom)
library(readr)
library(httr)
library(jsonlite)

### Load in GnomAD data.
########################
## GnomAD data was downloaded for the whole KCNH2 gene on the 26/04/2025. 
# Note, the downloaded data excludes failed variants, except those AC0. 
KCNH2_GnomAD <- read.csv("gnomAD_v4.1.1_ENSG00000055118_2026_04_26_14_38_34.csv")

### Load in the variant data.
#############################
Results_variant_recoder_df <- read.csv("Results_variant_recoder_protein_df_final.csv")

### Match variants with GnomAD data.
####################################
Results_variant_recoder_df2 <- Results_variant_recoder_df
Results_variant_recoder_df2$Final_cnomen <- gsub("NM_000238\\.4:", "", Results_variant_recoder_df2$Final_cnomen)

Results_variant_recoder_df2$Final_cnomen <- trimws(Results_variant_recoder_df2$Final_cnomen)
KCNH2_GnomAD$Transcript.Consequence <- trimws(KCNH2_GnomAD$Transcript.Consequence)

Results_variant_recoder_df2 <- Results_variant_recoder_df2 %>%
  mutate(Match_in_GnomAD = Final_cnomen %in% KCNH2_GnomAD$Transcript.Consequence)

Variants_merged_with_GnomAD <- merge(Results_variant_recoder_df2, KCNH2_GnomAD, by.x = "Final_cnomen", by.y = "Transcript.Consequence", all.x = TRUE)

Variant_list_no_GnomAD_data <- Results_variant_recoder_df2[Results_variant_recoder_df2$Match_in_GnomAD == FALSE,]
Variant_list_no_GnomAD_data <- unique(Variant_list_no_GnomAD_data$Final_cnomen)

## Loop through variants and determine classification based on GnomAD data. If 
## GnomAD data not available - inspect.

Allele_freq_criterion_result <- data.frame(HGVSC = character(), HGVSG = character(), HGVSG_alt = character(), HGVSP = character(), ID = character(),
                                           PM2_Supporting = logical(), BS1 = logical(), BA1 = logical(), AC_FAFpop = numeric(), AN_FAFpop = numeric(), 
                                           FAF_FAFpop = numeric(), LowerCI_FAFpop = numeric(), UpperCI_FAFpop = numeric(), Pop = character(), 
                                           Total_homozygotes = numeric(), Joint_filter = character(), Comment = character())

for (i in seq_len(nrow(Variants_merged_with_GnomAD))) {
  
  Variant_GnomAD_row <- Variants_merged_with_GnomAD[i, ] 
  
  HGVSC <- Variant_GnomAD_row$Final_cnomen
  HGVSG <- Variant_GnomAD_row$Final_gnomen
  HGVSG_alt <- Variant_GnomAD_row$Final_altgnomen
  HGVSP <- Variant_GnomAD_row$Final_pnomen
  Homozygote_counts <- Variant_GnomAD_row$Homozygote.Count 
  ID <- Variant_GnomAD_row$ID
  Joint_filter <- Variant_GnomAD_row$Filters...joint
  
  if (!is.na(Variant_GnomAD_row$gnomAD.ID)){
    
    # Ancestries, excluding those deemed as bottleneck according to GnomAD (Amish, Ashkenazi Jewish, Finnish, Remaining)
    Ancestries <- c("African.African.American", "Admixed.American", "East.Asian", "Middle.Eastern", "European..non.Finnish.", "South.Asian")
    Temp_FAF <- data.frame(HGVSC = character(), HGVSG = character(), HGVSG_alt = character(), HGVSP = character(), ID = character(), Pop = character(),
                           AC = numeric(), AN = numeric(), FAF = numeric(),stringsAsFactors = FALSE)
    
    for (Ancestry in Ancestries) {
      Column_name_AC <- paste0("Allele.Count.",Ancestry)
      Column_name_AN <- paste0("Allele.Number.",Ancestry)
      AC <- Variant_GnomAD_row[[Column_name_AC]]        # Allele Count.
      AN <- Variant_GnomAD_row[[Column_name_AN]]        # Allele Number.
      
      if (AC != 0) {
        FAF <- AC/AN                                    # Allele freq.
      } else if (AC == 0) {
        FAF <- NA
      }
      
      new_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID, Pop = Ancestry,
                            AC = AC,
                            AN = AN,
                            FAF = FAF)
      Temp_FAF <- rbind(Temp_FAF, new_row)
    } 
    
    # Need to capture the exome vs genome difs? 
    
    if (sum(Temp_FAF$AC) != 0) {
      highest_row <- Temp_FAF[which.max(Temp_FAF$FAF), ]
      Upper_CI <- binom.test(highest_row$AC, highest_row$AN, conf.level = 0.90)$conf.int[2]
      Lower_CI <- binom.test(highest_row$AC, highest_row$AN, conf.level = 0.90)$conf.int[1]
      
      PM2_threshold <- 0.00001  
      BS1_threshold <- 0.00005
      BA1_threshold <- 0.001
      
      if (Upper_CI <= PM2_threshold) {
        final_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID,
                                PM2_Supporting = TRUE, BS1 = FALSE, BA1 = FALSE, AC_FAFpop = highest_row$AC, AN_FAFpop = highest_row$AN, 
                                FAF_FAFpop = highest_row$FAF, LowerCI_FAFpop = Lower_CI, UpperCI_FAFpop = Upper_CI, Pop = highest_row$Pop, 
                                Total_homozygotes = Homozygote_counts, Joint_filter = Joint_filter, Comment = c("PM2_sup applied."))
        Allele_freq_criterion_result <- rbind(Allele_freq_criterion_result, final_row) 
      } else if (Lower_CI >= BA1_threshold) {
        final_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID,
                                PM2_Supporting = FALSE, BS1 = FALSE, BA1 = TRUE, AC_FAFpop = highest_row$AC, AN_FAFpop = highest_row$AN, 
                                FAF_FAFpop = highest_row$FAF, LowerCI_FAFpop = Lower_CI, UpperCI_FAFpop = Upper_CI, Pop = highest_row$Pop, 
                                Total_homozygotes = Homozygote_counts, Joint_filter = Joint_filter, Comment = c("BA1 applied."))
        Allele_freq_criterion_result <- rbind(Allele_freq_criterion_result, final_row) 
      } else if (Lower_CI >= BS1_threshold) {
        final_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID,
                                PM2_Supporting = FALSE, BS1 = TRUE, BA1 = FALSE, AC_FAFpop = highest_row$AC, AN_FAFpop = highest_row$AN, 
                                FAF_FAFpop = highest_row$FAF, LowerCI_FAFpop = Lower_CI, UpperCI_FAFpop = Upper_CI, Pop = highest_row$Pop, 
                                Total_homozygotes = Homozygote_counts, Joint_filter = Joint_filter, Comment = c("BS1 applied."))
        Allele_freq_criterion_result <- rbind(Allele_freq_criterion_result, final_row) 
      } else {
        final_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID,
                                PM2_Supporting = FALSE, BS1 = FALSE, BA1 = FALSE, AC_FAFpop = highest_row$AC, AN_FAFpop = highest_row$AN, 
                                FAF_FAFpop = highest_row$FAF, LowerCI_FAFpop = Lower_CI, UpperCI_FAFpop = Upper_CI, Pop = highest_row$Pop, 
                                Total_homozygotes = Homozygote_counts, Joint_filter = Joint_filter, Comment = c("Allele freq data fits no criteria."))
        Allele_freq_criterion_result <- rbind(Allele_freq_criterion_result, final_row) 
      }
    } else if (sum(Temp_FAF$AC) == 0) {
      final_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID,
                              PM2_Supporting = NA, BS1 = NA, BA1 = NA, AC_FAFpop = NA, AN_FAFpop = NA, 
                              FAF_FAFpop = NA, Pop = NA, LowerCI_FAFpop = NA, UpperCI_FAFpop = NA, 
                              Total_homozygotes = Homozygote_counts, Joint_filter = Joint_filter, Comment = c("No non-bottleneck pop with frequency data."))
      Allele_freq_criterion_result <- rbind(Allele_freq_criterion_result, final_row)
    }
  } else if (is.na(Variant_GnomAD_row$gnomAD.ID)) {
    final_row <- data.frame(HGVSC = HGVSC, HGVSG = HGVSG, HGVSG_alt = HGVSG_alt, HGVSP = HGVSP, ID = ID,
                            PM2_Supporting = NA, BS1 = NA, BA1 = NA, AC_FAFpop = NA, AN_FAFpop = NA, 
                            FAF_FAFpop = NA, Pop = NA, LowerCI_FAFpop = NA, UpperCI_FAFpop = NA, 
                            Total_homozygotes = Homozygote_counts, Joint_filter = NA, Comment = c("No pop with frequency data. PM2 query."))
    Allele_freq_criterion_result <- rbind(Allele_freq_criterion_result, final_row)
    message("Variant ", HGVSC, " has no available GnomAD data. To inspect.")
  }
}

## Two variants (ID 239, 574) were further examined on Alamut as they apparently have no non-bottleneck pop with frequency data available. 
## Both only have a single occurrence in the European Finnish population. Therefore, both meet PM2, assign this. 
Allele_freq_criterion_result$PM2_Supporting[Allele_freq_criterion_result$ID %in% c("239", "574")] <- TRUE
Allele_freq_criterion_result$BS1[Allele_freq_criterion_result$ID %in% c("239", "574")] <- FALSE
Allele_freq_criterion_result$BA1[Allele_freq_criterion_result$ID %in% c("239", "574")] <- FALSE

## For variants where there is no pop frequency data and a PM2 query, change PM2_Supporting to TRUE and BS1/BA1 to FALSE.
no_pop_rows <- Allele_freq_criterion_result$Comment == "No pop with frequency data. PM2 query."
Allele_freq_criterion_result$PM2_Supporting[no_pop_rows] <- TRUE
Allele_freq_criterion_result$BS1[no_pop_rows]            <- FALSE
Allele_freq_criterion_result$BA1[no_pop_rows]            <- FALSE

## For entries which met PM2 (and pop data exists in GnomAD) ensure there are no homozygotes in pop which could counter this evidence. 
PM2_homozygotes <- Allele_freq_criterion_result[Allele_freq_criterion_result$PM2_Supporting == TRUE &
                                                Allele_freq_criterion_result$Total_homozygotes > 0 &
                                                !is.na(Allele_freq_criterion_result$Total_homozygotes),]
# No cases where PM2 has been applied and there are homozygotes in the population. 

## Save this result. 
write.csv(Allele_freq_criterion_result, "Allele_freq_criterion_result_GnomAD_260426.csv", row.names = FALSE, quote = TRUE)
