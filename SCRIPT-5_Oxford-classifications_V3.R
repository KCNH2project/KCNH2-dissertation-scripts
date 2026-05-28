### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load in required libraries.
###############################
library(DBI)
library(RSQLite)
library(rvest)
library(readr) 
library(dplyr)
library(lubridate)

### Create file paths.
######################
## Variant list file path.
variant_file <- paste0(getwd(),"/All_variants_cnomen.csv")

## Alamut db path.
# File copied into C drive on the 22/11/2025. 
db_path <- paste0(getwd(),"/Variant database/Variant_Constitutional_AVPv1-12 10-03-25.db")

## Connect to database
con <- dbConnect(RSQLite::SQLite(), db_path)

## Check database structure to form rest of code. 
tables <- dbListTables(con)
print(tables)

columns <- dbListFields(con, "variant_patient")
print(columns)

## Read variants
variants <- read.csv(variant_file, header = FALSE, stringsAsFactors = FALSE)[,1]

# Change nomenclature of one variant temporarily, as otherwise won't search on Alamut. 
variants[variants == "c.982_983insATTATGCGCTACC"] <- "c.982_983ins13"

### Collect Alamut data for each variant in list.
#################################################
## Create empty data frame to store variant information from Alamut. 
Variant_OX_classifications <- data.frame(variant = "", occurence_ID = "", comment = "", comment_created = "", comment_createdby = "", comment_updated = "", comment_updatedby = "", 
                                         classification = "", classification_created = "", classification_updated = "", classification_updatedby = "")

## Loop over variants
for (i in seq_along(variants)) {
  variant <- variants[i]
  
  # Progress message for every variant
  message(paste("Processing variant", i, "/", length(variants), ":", variant))
  
  tryCatch({
    # Reset all variables at the start of the loop
    res_variantid <- data.frame()
    res_comment <- data.frame()
    res_class <- data.frame()
    variant_id <- NA
    classification <- ""
    classification_created_recent <- ""
    classification_updated_recent <- ""
    classification_updatedby_recent <- ""
    temp_df <- data.frame()
    
    # Get variant_id
    res_variantid <- dbGetQuery(con, 
                                'SELECT variant_id 
                                 FROM gene_annotation 
                                 WHERE gene_symbol = "KCNH2" 
                                 AND cNomen = ?', 
                                params = list(variant))
    
    if (nrow(res_variantid) == 0) {
      stop("Variant not found in gene_annotation")
    }
    variant_id <- res_variantid$variant_id[1]
    
    # Get comments + metadata in one query
    res_comment <- dbGetQuery(con, 
                              'SELECT patient_id, comment, created, updated, created_by, updated_by 
                               FROM variant_patient 
                               WHERE variant_id = ?', 
                              params = list(variant_id))
    
    if (nrow(res_comment) > 0) {
      # Clean HTML -> text
      res_comment$comment <- sapply(res_comment$comment, function(html_text) {
        text <- xml2::read_html(paste0("<body>", html_text, "</body>")) |> rvest::html_text()
        text <- gsub("\n", " | ", text)
        sub("^ \\|  \\| p, li \\{ white-space: pre-wrap; \\} \\| ", "", text)
      })
    }
    
    # Get classification history
    res_class <- dbGetQuery(con, 
                            'SELECT classification, created, updated, updated_by 
                             FROM variant_history 
                             WHERE variant_id = ?', 
                            params = list(variant_id))
    
    classification <- if (nrow(res_class) > 0) {
      tail(res_class$classification, 1)
    } else {
      ""
    }
    
    classification_created_recent <- if (nrow(res_class) > 0) {
      tail(res_class$created, 1)
    } else {
      ""
    }
    
    classification_updated_recent <- if (nrow(res_class) > 0) {
      tail(res_class$updated, 1)
    } else {
      ""
    }
    
    classification_updatedby_recent <- if (nrow(res_class) > 0) {
      tail(res_class$updated_by, 1)
    } else {
      ""
    }
    
    # Write rows
    if (nrow(res_comment) > 0) {
      for (j in seq_len(nrow(res_comment))) {
        temp_df <- data.frame(
          variant = variant,
          occurence_ID = res_comment$patient_id[j],
          comment = res_comment$comment[j],
          comment_created = res_comment$created[j],
          comment_createdby = res_comment$created_by[j],
          comment_updated = res_comment$updated[j],
          comment_updatedby = res_comment$updated_by[j],
          classification = classification,
          classification_created = classification_created_recent,
          classification_updated = classification_updated_recent,
          classification_updatedby = classification_updatedby_recent
        )
        Variant_OX_classifications <- rbind(Variant_OX_classifications, temp_df)
      }
    } else {
      temp_df <- data.frame(
        variant = variant,
        occurence_ID = "",
        comment = "",
        comment_created = "",
        comment_createdby = "",
        comment_updated = "",
        comment_updatedby = "",
        classification = classification,
        classification_created = classification_created_recent,
        classification_updated = classification_updated_recent,
        classification_updatedby = classification_updatedby_recent
      )
      Variant_OX_classifications <- rbind(Variant_OX_classifications, temp_df)
    }
    
  }, error = function(e) {
    message(paste("Error for variant", variant, ":", e$message))
    temp_df <- data.frame(
      variant = variant,
      occurence_ID = "",
      comment = "Not recorded in Oxford Alamut database.",
      comment_created = "",
      comment_createdby = "",
      comment_updated = "",
      comment_updatedby = "",
      classification = "",
      classification_created = "",
      classification_updated = "",
      classification_updatedby = ""
    )
    Variant_OX_classifications <<- rbind(Variant_OX_classifications, temp_df)
  })
}

# Close DB
dbDisconnect(con)

### Check all Oxford variants provided have Alamut data. 
########################################################
## Oxford variants previously established.
Results_variant_recoder_df <- read.csv("Results_variant_recoder_protein_df_final.csv")
Oxford_variant_recoder_final <- Results_variant_recoder_df[Results_variant_recoder_df$Lab == "Oxford",]
Oxford_variant_recoder_final$Final_cnomen <- gsub("NM_000238\\.4:", "", Oxford_variant_recoder_final$Final_cnomen)
Oxford_variant_list <- Oxford_variant_recoder_final$Final_cnomen

## Compare to the classifications successfully pulled from the Oxford Alamut database. 
Variant_OX_classifications_success <- Variant_OX_classifications[!Variant_OX_classifications$comment == "Not recorded in Oxford Alamut database.",]
# Change nomenclature of one variant again, as otherwise won't be found in cross-comparison. 
Variant_OX_classifications_success$variant[Variant_OX_classifications_success$variant == "c.982_983ins13"] <- "c.982_983insATTATGCGCTACC"
Oxford_variant_list_df <- data.frame(Variants = Oxford_variant_list, In_Alamut = "")
Oxford_variant_list_df$In_Alamut <- Oxford_variant_list_df$Variants %in% Variant_OX_classifications_success$variant

### Check Oxford KCNH2 variant list is up-to-date. 
##################################################
## Connect to database
con <- dbConnect(RSQLite::SQLite(), db_path)

## Retrieve all KCNH2 variants in Oxford database. 
All_Ox_KCNH2_var <- dbGetQuery(con, 
                   'SELECT cNomen 
                    FROM gene_annotation 
                    WHERE gene_symbol = "KCNH2"')

## Cross compare to existing Oxford KCNH2 variant list, extract those which are new. 
All_Ox_KCNH2_var_new <- All_Ox_KCNH2_var[!All_Ox_KCNH2_var$cNomen %in% Oxford_variant_list_df$Variants,]
# Remove the c.982_983ins13 variant, as this isn't new it's just differences in naming. 
All_Ox_KCNH2_var_new <- All_Ox_KCNH2_var_new[!All_Ox_KCNH2_var_new == "c.982_983ins13"]

## Create empty df to hold classifications for new Oxford variants. 
All_Ox_KCNH2_var_classifications <- data.frame(variant = "", occurence_ID = "", comment = "", comment_created = "", comment_createdby = "", comment_updated = "", comment_updatedby = "", 
                                               classification = "", classification_created = "", classification_updated = "", classification_updatedby = "")

## Loop over variants
for (i in seq_along(All_Ox_KCNH2_var_new)) {
  Ox_variant <- All_Ox_KCNH2_var_new[i]
  
  # Progress message for every variant
  message(paste("Processing variant", i, "/", length(All_Ox_KCNH2_var_new), ":", Ox_variant))
  
  tryCatch({
    # Reset all variables at the start of the loop
    res_variantid <- data.frame()
    res_comment <- data.frame()
    res_class <- data.frame()
    variant_id <- NA
    classification <- ""
    classification_created_recent <- ""
    classification_updated_recent <- ""
    classification_updatedby_recent <- ""
    temp_df <- data.frame()
    
    # Get variant_id
    res_variantid <- dbGetQuery(con, 
                                'SELECT variant_id 
                                 FROM gene_annotation 
                                 WHERE gene_symbol = "KCNH2" 
                                 AND cNomen = ?', 
                                params = list(Ox_variant))
    
    if (nrow(res_variantid) == 0) {
      stop("Variant not found in gene_annotation")
    }
    variant_id <- res_variantid$variant_id[1]
    
    # Get comments + metadata in one query
    res_comment <- dbGetQuery(con, 
                              'SELECT patient_id, comment, created, updated, created_by, updated_by 
                               FROM variant_patient 
                               WHERE variant_id = ?', 
                              params = list(variant_id))
    
    if (nrow(res_comment) > 0) {
      # Clean HTML -> text
      res_comment$comment <- sapply(res_comment$comment, function(html_text) {
        text <- xml2::read_html(paste0("<body>", html_text, "</body>")) |> rvest::html_text()
        text <- gsub("\n", " | ", text)
        sub("^ \\|  \\| p, li \\{ white-space: pre-wrap; \\} \\| ", "", text)
      })
    }
    
    # Get classification history
    res_class <- dbGetQuery(con, 
                            'SELECT classification, created, updated, updated_by 
                             FROM variant_history 
                             WHERE variant_id = ?', 
                            params = list(variant_id))
    
    classification <- if (nrow(res_class) > 0) {
      tail(res_class$classification, 1)
    } else {
      ""
    }
    
    classification_created_recent <- if (nrow(res_class) > 0) {
      tail(res_class$created, 1)
    } else {
      ""
    }
    
    classification_updated_recent <- if (nrow(res_class) > 0) {
      tail(res_class$updated, 1)
    } else {
      ""
    }
    
    classification_updatedby_recent <- if (nrow(res_class) > 0) {
      tail(res_class$updated_by, 1)
    } else {
      ""
    }
    
    # Write rows
    if (nrow(res_comment) > 0) {
      for (j in seq_len(nrow(res_comment))) {
        temp_df <- data.frame(
          variant = Ox_variant,
          occurence_ID = res_comment$patient_id[j],
          comment = res_comment$comment[j],
          comment_created = res_comment$created[j],
          comment_createdby = res_comment$created_by[j],
          comment_updated = res_comment$updated[j],
          comment_updatedby = res_comment$updated_by[j],
          classification = classification,
          classification_created = classification_created_recent,
          classification_updated = classification_updated_recent,
          classification_updatedby = classification_updatedby_recent
        )
        All_Ox_KCNH2_var_classifications <- rbind(All_Ox_KCNH2_var_classifications, temp_df)
      }
    } else {
      temp_df <- data.frame(
        variant = Ox_variant,
        occurence_ID = "",
        comment = "",
        comment_created = "",
        comment_createdby = "",
        comment_updated = "",
        comment_updatedby = "",
        classification = classification,
        classification_created = classification_created_recent,
        classification_updated = classification_updated_recent,
        classification_updatedby = classification_updatedby_recent
      )
      All_Ox_KCNH2_var_classifications <- rbind(All_Ox_KCNH2_var_classifications, temp_df)
    }
    
  }, error = function(e) {
    message(paste("Error for variant", Ox_variant, ":", e$message))
    temp_df <- data.frame(
      variant = Ox_variant,
      occurence_ID = "",
      comment = "Not recorded in Oxford Alamut database.",
      comment_created = "",
      comment_createdby = "",
      comment_updated = "",
      comment_updatedby = "",
      classification = "",
      classification_created = "",
      classification_updated = "",
      classification_updatedby = ""
    )
    All_Ox_KCNH2_var_classifications <<- rbind(All_Ox_KCNH2_var_classifications, temp_df)
  })
}

# Close DB
dbDisconnect(con)

### Gather classification and criteria for each variant (if possible). 
######################################################################
## Fix date and time columns. 
# Parse as POSIXct
comment_created_parsed <- ymd_hms(Variant_OX_classifications_success$comment_created)
comment_updated_parsed <- ymd_hms(Variant_OX_classifications_success$comment_updated)
classification_created_parsed <- ymd_hms(Variant_OX_classifications_success$classification_created)
classification_updated_parsed <- ymd_hms(Variant_OX_classifications_success$classification_updated)

# Combine into a data.frame
Variant_OX_classifications_success_datetimefix <- Variant_OX_classifications_success
Variant_OX_classifications_success_datetimefix$comment_created <- comment_created_parsed
Variant_OX_classifications_success_datetimefix$comment_updated <- comment_updated_parsed
Variant_OX_classifications_success_datetimefix$classification_created <- classification_created_parsed
Variant_OX_classifications_success_datetimefix$classification_updated <- classification_updated_parsed

typeof(Variant_OX_classifications_success_datetimefix$comment_created)

## Create a data frame which will need to be manually edited to select the most recent criterion and classification
## decided for each variant in the list in Oxford. 
Oxford_classification_all_variants <- data.frame(Variant_cnomen = "", 
                                                 PVS1 = "", PS1 = "", PS2 = "", PS3 = "", PS4 = "", PM1 = "", PM4 = "", PM5 = "", PM6 = "", 
                                                 PP1 = "", PP2 = "", PP3 = "", PM2 = "", BA1 = "", BS1 = "", BS3 = "", BS4 = "", BP2 = "", 
                                                 BP4 = "", BP7 = "", Final_classification = "", Own_comment = "", Classified_date = "", Alamut_classification = "", 
                                                 Alamut_classification_created = as.POSIXct("2018-05-26 11:24:52", origin = "1970-01-01", tz = "UTC"), 
                                                 Alamut_classification_updated = as.POSIXct("2018-05-26 11:24:52", origin = "1970-01-01", tz = "UTC"), 
                                                 Alamut_classification_updatedby = "")

# For loop to go through each variant and input the variant cnomen & classification details.
for (cnomen1 in unique(Variant_OX_classifications_success_datetimefix$variant)) {
  temp_df1 <- Variant_OX_classifications_success_datetimefix[Variant_OX_classifications_success_datetimefix$variant == cnomen1,]
  Oxford_classification_all_variants_row <- data.frame(Variant_cnomen = cnomen1, 
                                                   PVS1 = "", PS1 = "", PS2 = "", PS3 = "", PS4 = "", PM1 = "", PM4 = "", PM5 = "", PM6 = "", 
                                                   PP1 = "", PP2 = "", PP3 = "", PM2 = "", BA1 = "", BS1 = "", BS3 = "", BS4 = "", BP2 = "", 
                                                   BP4 = "", BP7 = "", 
                                                   Final_classification = "", Own_comment = "", Classified_date = "",
                                                   Alamut_classification = unique(temp_df1$classification), 
                                                   Alamut_classification_created = unique(temp_df1$classification_created), 
                                                   Alamut_classification_updated = unique(temp_df1$classification_updated), 
                                                   Alamut_classification_updatedby = unique(temp_df1$classification_updatedby))
  Oxford_classification_all_variants <- rbind(Oxford_classification_all_variants, Oxford_classification_all_variants_row)
}

# Remove empty rows from df. 
Oxford_classification_all_variants <- Oxford_classification_all_variants[-c(1,2),]

## Save this df for manual editing. 
write.csv(Oxford_classification_all_variants, "Oxford_classifications_for_all_lab_variants.csv", row.names = FALSE, quote = TRUE)

## Select the most recent top 5 rows based on the comment updated date/time. 
# Create empty df with same structure as Variant_OX_classifications_success_datetimefix.
Variant_OX_classifications_top5_condensed <- data.frame(variant = "", occurence_ID = "", comment = "", comment_created = "", comment_createdby = "", comment_updated = "", 
                                                        comment_updatedby = "", classification = "", classification_created = "", classification_updated = "", 
                                                        classification_updatedby = "")

# For loop to select 5 most recent occurrences for each variant.
for (cnomen in unique(Variant_OX_classifications_success_datetimefix$variant)) {
  Variant_alamut_info <- data.frame()
  Variant_alamut_info <- Variant_OX_classifications_success_datetimefix[Variant_OX_classifications_success_datetimefix$variant == cnomen,]
  Variant_alamut_info <- Variant_alamut_info[order(Variant_alamut_info$comment_updated, decreasing = TRUE),]
  Variant_alamut_info <- head(Variant_alamut_info, 5)
  Variant_alamut_info <- Variant_alamut_info %>%
    group_by(variant) %>%
    summarise(
      across(
        everything(),
        ~ {
          vals <- trimws(na.omit(.x))              # remove NA + strip whitespace
          vals <- vals[vals != ""]                 # drop empty strings
          if (length(vals) == 0) NA_character_ else paste(unique(vals), collapse = " END OF ENTRY \n")
        }
      ),
      .groups = "drop"
    )
  Variant_OX_classifications_top5_condensed <- rbind(Variant_OX_classifications_top5_condensed, Variant_alamut_info)
}

# Remove empty rows from df. 
Variant_OX_classifications_top5_condensed <- Variant_OX_classifications_top5_condensed[-c(1,2),]

## Save this df for manual editing. Final classification determined using 2024 guidelines. 
write.csv(Variant_OX_classifications_top5_condensed, "All_variant_OX_classifications_top5_condensed.csv", row.names = FALSE, quote = TRUE)
