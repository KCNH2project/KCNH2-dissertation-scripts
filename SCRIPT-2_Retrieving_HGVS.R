### Set working directory.
##########################
setwd("C:/Users/chloe.greve/OneDrive - Oxford University Hospitals NHS Foundation Trust/Documents/Research project/Project/Variants/R work/Updated final final scripts")

### Load in required libraries.
###############################
library(tidyr)
library(dplyr)
library(stringr)
library(httr)
library(jsonlite)
library(xml2)

### Load in the variants given by the labs.
###########################################
Collated_variants <- read.csv("Collated_variant_list.csv", row.names = NULL)

## Add an ID column to keep rows unique when processing. 
Collated_variants$ID <- seq_len(nrow(Collated_variants))

## Remove white space around the c. and p. nomenclature.
Collated_variants$c.nomenclature <- str_squish(Collated_variants$c.nomenclature)
Collated_variants$p.nomenclature <- str_squish(Collated_variants$p.nomenclature)
Collated_variants$g.nomenclature <- str_squish(Collated_variants$g.nomenclature)

## Add up-to-date mRNA transcript and corresponding protein reference sequence. 
Collated_variants$c.nomenclature <- paste0("NM_000238.4:", Collated_variants$c.nomenclature)
Collated_variants$p.nomenclature <- paste0("NP_000229.1:", Collated_variants$p.nomenclature)


### Determine number of variants working with prior to filtering. 
#################################################################

## Whole data set.
length(unique(Collated_variants$c.nomenclature))

Total_missense_variants <- Collated_variants[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Collated_variants$p.nomenclature
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Collated_variants$p.nomenclature
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Collated_variants$c.nomenclature
    ),
]

length(unique(Total_missense_variants$c.nomenclature))

## Aberdeen.
Aberdeen_df <- Collated_variants[Collated_variants$Lab == "Aberdeen", ]

# Determine number of unique variants. 
length(unique(Aberdeen_df$c.nomenclature))

Aberdeen_missense_variants <- Aberdeen_df[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Aberdeen_df$p.nomenclature
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Aberdeen_df$p.nomenclature
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Aberdeen_df$c.nomenclature
    ),
]

length(unique(Aberdeen_missense_variants$c.nomenclature))


## Belfast.
Belfast_df <- Collated_variants[Collated_variants$Lab == "Belfast", ]

# Determine number of unique variants. 
length(unique(Belfast_df$c.nomenclature))

Belfast_missense_variants <- Belfast_df[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Belfast_df$p.nomenclature
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Belfast_df$p.nomenclature
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Belfast_df$c.nomenclature
    ),
]

length(unique(Belfast_missense_variants$c.nomenclature))


## Brompton.
Brompton_df <- Collated_variants[Collated_variants$Lab == "Brompton", ]

# Determine number of unique variants. 
length(unique(Brompton_df$c.nomenclature))

Brompton_missense_variants <- Brompton_df[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Brompton_df$p.nomenclature
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Brompton_df$p.nomenclature
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Brompton_df$c.nomenclature
    ),
]

length(unique(Brompton_missense_variants$c.nomenclature))


## Manchester.
Manchester_df <- Collated_variants[Collated_variants$Lab == "Manchester", ]

# Determine number of unique variants. 
length(unique(Manchester_df$c.nomenclature))

Manchester_missense_variants <- Manchester_df[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Manchester_df$p.nomenclature
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Manchester_df$p.nomenclature
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Manchester_df$c.nomenclature
    ),
]

length(unique(Manchester_missense_variants$c.nomenclature))


## Oxford.
Oxford_df <- Collated_variants[Collated_variants$Lab == "Oxford", ]

# Determine number of unique variants. 
length(unique(Oxford_df$c.nomenclature))

Oxford_missense_variants <- Oxford_df[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Oxford_df$p.nomenclature
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Oxford_df$p.nomenclature
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Oxford_df$c.nomenclature
    ),
]

length(unique(Oxford_missense_variants$c.nomenclature))

### Set up structure to send queries to Variant Recoder using c. and p.
#######################################################################

## Define the server (i.e. where variant recoder is held).
server <- "https://rest.ensembl.org"

## Create empty list to hold results using c. search.  
Results_variant_recoder_list <- list()

## Define retry parameters
max_retries <- 5  # Maximum number of retries
retry_delay <- 5   # Delay between retries in seconds

## Send variant queries to Variant Recoder. 
for (i in seq_len(nrow(Collated_variants))) {
  # Go through each row/variant and retrieve VEP info. 
  Collated_variants_row <- Collated_variants[i, ] 
  
  # Initialize retry counter
  attempt <- 1
  success <- FALSE
  skip_chunk <- FALSE
  
  Variant <- Collated_variants_row$c.nomenclature
  print(paste("Attempting", Variant ,"HGVS retrieval."))
  
  while (attempt <= max_retries && !success && !skip_chunk) { 
    tryCatch({
      # Make the GET request
      ext <- paste0("/variant_recoder/human/",Variant,"?vcf_string=1")
      r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
      
      # Error check. 
      stop_for_status(r) 
      
      # Parse and store response. 
      response_content <- content(r)
      
      # Collect HGVS g., c., and p. for each variant. 
      if (grepl("c\\.[+-]?\\d+(?:[+-]\\d+)?_[+-]?\\d+(?:[+-]\\d+)?del", Variant)) { 
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        
        if (!is.na(matches[1])) {
          alt_position <- as.numeric(matches[2])
          alt_position_calc <- as.numeric(alt_position) + 1
          sequence <- matches[3]
          alt_sequence_calc <- nchar(sequence) - 1 + alt_position_calc
        } else if (is.na(matches[1])) {
          hgvsg_alt_pre <- response_content[[1]][[Link]]$vcf_string[[1]]
          matches <- str_match(hgvsg_alt_pre, "^([0-9XYM]+)-([0-9]+)-([ACGT]+)-([ACGT]+)$")
          alt_position <- as.numeric(matches[3])
          alt_position_calc <- as.numeric(alt_position) + 1
          sequence <- matches[4]
          alt_sequence_calc <- nchar(sequence) - 2 + alt_position_calc
        }
        
        g_first_number <- as.numeric(sub(".*g\\.(\\d+)_.*", "\\1", hgvsg))
        
        if (as.numeric(g_first_number) >= (alt_position + 2)) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_calc), 
                              "_", paste0(alt_sequence_calc), "del") ## If delins need to add nucleotide count after delins. 
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+_\\d+delins", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+):{1,2}([ACGT]+)")
        alt_position <- as.numeric(matches[2])
        alt_position_calc <- as.numeric(alt_position) + 1
        sequence_del <- matches[3]
        sequence_ins <- matches[4]
        alt_sequence_calc <- nchar(sequence_del) - 1 + alt_position_calc
        g_first_number <- as.numeric(sub(".*g\\.(\\d+)_.*", "\\1", hgvsg))
        
        if (as.numeric(g_first_number) >= (alt_position + 2)) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_calc), 
                              "_", paste0(alt_sequence_calc), "delins", sequence_ins)  
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\.\\d+(?:[+-]\\d+)?del(?![_\\d])", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position <- as.numeric(matches[2])
        alt_position_calc <- as.numeric(alt_position) + 1
        g_first_number <- as.numeric(sub(".*g\\.(\\d+).*", "\\1", hgvsg))
        
        if (as.numeric(g_first_number) >= (alt_position + 2)) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_calc),"del")
        } else {
          hgvsg_alt <-"No alt."
        }
        
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+_\\d+dup", Variant)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position_end <- as.numeric(matches[2])
        sequence <- matches[3]
        alt_position_start <- alt_position_end - (nchar(sequence) - 1)
        
        # For dups >1 bp assuming always alt as repetitive region - to assess when doing variant 
        # interpretation. 
        hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_start), 
                            "_", paste0(alt_position_end), "dup")
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+dup(?![_\\d])", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position <- as.numeric(matches[2])
        g_first_number <- as.numeric(sub(".*g\\.(\\d+).*", "\\1", hgvsg))
        
        if (g_first_number != alt_position) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position),"dup")
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\.[\\d\\-]+_[\\d\\-]+ins(?![a-z]*delins)", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position_start <- as.numeric(matches[2])
        alt_position_end <- alt_position_start + 1
        g_first_number <- as.numeric(sub(".*g\\.(\\d+)_.*", "\\1", hgvsg))
        ins <- matches[3]
        
        if (g_first_number != alt_position_start) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_start), 
                              "_", paste0(alt_position_end), "ins",ins) ## Edit this to include the nucleotide letter count.
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+[ACGT]>[ACGT]", Variant)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        hgvsg_alt <-"No alt."
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      }
      
      Results_variant_recoder_list[[Collated_variants_row$ID]] <-  data.frame(hgvsg = hgvsg, hgvsg_alt = hgvsg_alt, hgvsc = hgvsc, 
                                                                       hgvsp = hgvsp, ID = Collated_variants_row$ID)
      
      success <- TRUE  # Mark the request as successful
      
      cat("Processed entry",Variant, "\n")
    }, error = function(e) {
      # If an error occurs, print error message. 
      message("Error processing ", i, ": ", e$message)
      
      # If HTTP 400 error skip, as this means the request is incorrect. 
      if (grepl("Bad Request \\(HTTP 400\\)", e$message)) {
        skip_chunk <<- TRUE
        success <- FALSE
      }
      
      # If max retries have not been reached, wait and try again.
      if (attempt < max_retries && !skip_chunk) {
        message("Retrying... Attempt ", attempt + 1, " of ", max_retries)
        Sys.sleep(retry_delay)  # Wait before retrying
      }
      
      # If max retries are reached, fail the call. 
      if (attempt == max_retries) {
        success <- FALSE
      }
    })
    # Increment the retry attempt counter.
    attempt <- attempt + 1
  }
  
  # If all retry attempts fail, log in results list and continue.
  if (success == FALSE && !skip_chunk) {
    message("Failed to process entry ", i, " after ", max_retries, " attempts.")
    Results_variant_recoder_list[[Collated_variants_row$ID]] <- c("Unable to process variant - to inspect.", Collated_variants_row$c.nomenclature)
  }
  if (success == FALSE && skip_chunk) {
    message("Failed to process entry due to error (HTTP 400). Skipping.")
    Results_variant_recoder_list[[Collated_variants_row$ID]] <- c("Unable to process variant - to inspect.", Collated_variants_row$c.nomenclature)
  }
  
  # Pause between chunks to avoid throttling.
  Sys.sleep(5)
}

## Create large data frame containing all the gathered variant info. 
# First, filter out any non-data frames or malformed entries using function. 
Results_variant_recoder_list_tidy <- Filter(function(x) {
  is.data.frame(x) && ncol(x) > 0 && nrow(x) > 0
}, Results_variant_recoder_list)

# Now merge variants into one large data frame
Results_variant_recoder_df <- do.call(rbind, Results_variant_recoder_list_tidy)

# Save these results.
#write.csv(Results_variant_recoder_df, "Results_variant_recoder_df.csv", row.names = FALSE, quote = TRUE)

# Check reads in correctly. 
Results_variant_recoder_df <- read.csv("Results_variant_recoder_df.csv")

## Create empty list to hold results using p. search.  
Results_variant_recoder_list_protein <- list()

## Send protein queries to Variant Recoder. 
for (i in seq_len(nrow(Collated_variants))) {
  # Go through each row/variant and retrieve VEP info. 
  Collated_variants_row <- Collated_variants[i, ] 
  
  # Initialize retry counter
  attempt <- 1
  success <- FALSE
  skip_chunk <- FALSE
  
  Variant <- Collated_variants_row$c.nomenclature
  Protein <- Collated_variants_row$p.nomenclature
  print(paste("Attempting", Protein ,"HGVS retrieval."))
  
  while (attempt <= max_retries && !success && !skip_chunk) { 
    tryCatch({
      # Make the GET request
      ext <- paste0("/variant_recoder/human/",Protein,"?vcf_string=1")
      r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
      
      # Error check. 
      stop_for_status(r) 
      
      # Parse and store response. 
      response_content <- content(r)
      
      # Collect HGVS g., c., and p. for each variant. 
      if (grepl("c\\.[+-]?\\d+(?:[+-]\\d+)?_[+-]?\\d+(?:[+-]\\d+)?del", Variant)) { 
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        
        if (!is.na(matches[1])) {
          alt_position <- as.numeric(matches[2])
          alt_position_calc <- as.numeric(alt_position) + 1
          sequence <- matches[3]
          alt_sequence_calc <- nchar(sequence) - 1 + alt_position_calc
        } else if (is.na(matches[1])) {
          hgvsg_alt_pre <- response_content[[1]][[Link]]$vcf_string[[1]]
          matches <- str_match(hgvsg_alt_pre, "^([0-9XYM]+)-([0-9]+)-([ACGT]+)-([ACGT]+)$")
          alt_position <- as.numeric(matches[3])
          alt_position_calc <- as.numeric(alt_position) + 1
          sequence <- matches[4]
          alt_sequence_calc <- nchar(sequence) - 2 + alt_position_calc
        }
        
        g_first_number <- as.numeric(sub(".*g\\.(\\d+)_.*", "\\1", hgvsg))
        
        if (as.numeric(g_first_number) >= (alt_position + 2)) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_calc), 
                              "_", paste0(alt_sequence_calc), "del") ## If delins need to add nucleotide count after delins. 
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+_\\d+delins", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+):{1,2}([ACGT]+)")
        alt_position <- as.numeric(matches[2])
        alt_position_calc <- as.numeric(alt_position) + 1
        sequence_del <- matches[3]
        sequence_ins <- matches[4]
        alt_sequence_calc <- nchar(sequence_del) - 1 + alt_position_calc
        g_first_number <- as.numeric(sub(".*g\\.(\\d+)_.*", "\\1", hgvsg))
        
        if (as.numeric(g_first_number) >= (alt_position + 2)) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_calc), 
                              "_", paste0(alt_sequence_calc), "delins", sequence_ins)  
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\.\\d+(?:[+-]\\d+)?del(?![_\\d])", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position <- as.numeric(matches[2])
        alt_position_calc <- as.numeric(alt_position) + 1
        g_first_number <- as.numeric(sub(".*g\\.(\\d+).*", "\\1", hgvsg))
        
        if (as.numeric(g_first_number) >= (alt_position + 2)) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_calc),"del")
        } else {
          hgvsg_alt <-"No alt."
        }
        
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+_\\d+dup", Variant)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position_end <- as.numeric(matches[2])
        sequence <- matches[3]
        alt_position_start <- alt_position_end - (nchar(sequence) - 1)
        
        # For dups >1 bp assuming always alt as repetitive region - to assess when doing variant 
        # interpretation. 
        hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_start), 
                            "_", paste0(alt_position_end), "dup")
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+dup(?![_\\d])", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position <- as.numeric(matches[2])
        g_first_number <- as.numeric(sub(".*g\\.(\\d+).*", "\\1", hgvsg))
        
        if (g_first_number != alt_position) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position),"dup")
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\.[\\d\\-]+_[\\d\\-]+ins(?![a-z]*delins)", Variant, perl = TRUE)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        
        hgvsg_alt_pre <- response_content[[1]][[Link]]$spdi[[1]]
        matches <- str_match(hgvsg_alt_pre, "NC_\\d+\\.\\d+:(\\d+):{1,2}([ACGT]+)")
        alt_position_start <- as.numeric(matches[2])
        alt_position_end <- alt_position_start + 1
        g_first_number <- as.numeric(sub(".*g\\.(\\d+)_.*", "\\1", hgvsg))
        ins <- matches[3]
        
        if (g_first_number != alt_position_start) {
          hgvsg_alt <- paste0("NC_000007.14:g.", paste0(alt_position_start), 
                              "_", paste0(alt_position_end), "ins",ins) ## Edit this to include the nucleotide letter count.
        } else {
          hgvsg_alt <-"No alt."
        }
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      } else if (grepl("\\d+[ACGT]>[ACGT]", Variant)) {
        
        Link <- names(response_content[[1]])[1]
        hgvsg <- response_content[[1]][[Link]]$hgvsg[[1]]
        hgvsg_alt <-"No alt."
        
        hgvsp_pre <- response_content[[1]][[Link]]$hgvsp
        hgvsp <- grep("^NP_000229\\.1:p\\.", hgvsp_pre, value = TRUE)
        
        if (length(hgvsp) == 0) {
          hgvsp <- "None retrieved"
        }
        
        hgvsc_pre <- response_content[[1]][[Link]]$hgvsc
        hgvsc <- grep("^NM_000238\\.4:c\\.", hgvsc_pre, value = TRUE)
        
      }
      
      Results_variant_recoder_list_protein[[Collated_variants_row$ID]] <-  data.frame(hgvsg = hgvsg, hgvsg_alt = hgvsg_alt, hgvsc = hgvsc, 
                                                                              hgvsp = hgvsp, ID = Collated_variants_row$ID)
      
      success <- TRUE  # Mark the request as successful
      
      cat("Processed entry",Protein, "\n")
    }, error = function(e) {
      # If an error occurs, print error message. 
      message("Error processing ", i, ": ", e$message)
      
      # If HTTP 400 error skip, as this means the request is incorrect. 
      if (grepl("Bad Request \\(HTTP 400\\)", e$message)) {
        skip_chunk <- TRUE
        success <- FALSE
      }
      
      # If max retries have not been reached, wait and try again.
      if (attempt < max_retries && !skip_chunk) {
        message("Retrying... Attempt ", attempt + 1, " of ", max_retries)
        Sys.sleep(retry_delay)  # Wait before retrying
      }
      
      # If max retries are reached, fail the call. 
      if (attempt == max_retries) {
        success <- FALSE
      }
    })
    # Increment the retry attempt counter.
    attempt <- attempt + 1
  }
  
  # If all retry attempts fail, log in results list and continue.
  if (success == FALSE && !skip_chunk) {
    message("Failed to process entry ", i, " after ", max_retries, " attempts.")
    Results_variant_recoder_list_protein[[Collated_variants_row$ID]] <- c("Unable to process variant - to inspect.", Collated_variants_row$p.nomenclature)
  }
  if (success == FALSE && skip_chunk) {
    message("Failed to process entry due to error (HTTP 400). Skipping.")
    Results_variant_recoder_list_protein[[Collated_variants_row$ID]] <- c("Unable to process variant - to inspect.", Collated_variants_row$p.nomenclature)
  }
  
  # Pause between chunks to avoid throttling.
  Sys.sleep(5)
}

## Create large data frame containing all the gathered variant info. 
# First, filter out any non-data frames or malformed entries using function. 
Results_variant_recoder_list_protein_tidy <- Filter(function(x) {
  is.data.frame(x) && ncol(x) > 0 && nrow(x) > 0
}, Results_variant_recoder_list_protein)

# Now merge variants into one large data frame
Results_variant_recoder_protein_df <- do.call(rbind, Results_variant_recoder_list_protein_tidy)

# Change column names so can compare with c. nomenclature query. 
colnames(Results_variant_recoder_protein_df) <- c("hgvsg_protein", "hgvsg_alt_protein", "hgvsc_protein", 
                                                  "hgvsp_protein", "ID")

# Save these results.
#write.csv(Results_variant_recoder_protein_df, "Results_variant_recoder_protein_df.csv", row.names = FALSE, quote = TRUE)

# Check reads in correctly. 
Results_variant_recoder_protein_df <- read.csv("Results_variant_recoder_protein_df.csv")


### Compare results with provided nomenclature.
###############################################

## Merge the HGVS retrieved and the existing collated variant data frame by row ID. 
# Some of the protein queries pulled two versions of the c., remove duplicates. 
Results_variant_recoder_protein_df_unique <- Results_variant_recoder_protein_df %>%
                                             distinct(ID, .keep_all = TRUE)

merged_df <- left_join(Collated_variants, Results_variant_recoder_df, by = "ID") %>%
             left_join(Results_variant_recoder_protein_df_unique, by = "ID")

## Change lab-provided amino acid names to 3 letter names to match variant recoder results. 
aa_map <- c(A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys",
            E = "Glu", Q = "Gln", G = "Gly", H = "His", I = "Ile",
            L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro",
            S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val")

expand_variant <- function(variant) {
  # Skip if already has 3-letter codes
  if (grepl("p\\.[A-Z][a-z]{2}", variant)) {
    return(variant)
  }
  
  pattern <- "(?<=p\\.|_|[0-9])([A-Z])(?![a-z])"
  
  m <- gregexpr(pattern, variant, perl = TRUE)
  matches <- regmatches(variant, m)
  
  # Replace using aa_map
  matches <- lapply(matches, function(x) aa_map[x])
  
  regmatches(variant, m) <- matches
  variant
}

merged_df$p.nomenclature <- sapply(merged_df$p.nomenclature, expand_variant)

## Remove unnecessary nucleotide(s) at the end of duplication and deletion variants. 
merged_df$c.nomenclature <- sub("(dup)[ACGT]+$", "\\1", merged_df$c.nomenclature)
merged_df$c.nomenclature <- sub("(del)[ACGT]+$", "\\1", merged_df$c.nomenclature)

## Swap "*" for "Ter" in protein nomenclature.
merged_df$p.nomenclature <- sub("\\*", "Ter", merged_df$p.nomenclature)

## Conduct comparison of nomenclature provided vs retrieved for c. and p. using a T/F column.  
merged_df <- merged_df %>%
  mutate(Matching_nomenclature_coding_cquery = str_squish(merged_df$c.nomenclature) == str_squish(merged_df$hgvsc))

merged_df <- merged_df %>%
  mutate(Matching_nomenclature_protein_cquery  = str_squish(merged_df$p.nomenclature) == str_squish(merged_df$hgvsp))

merged_df <- merged_df %>%
  mutate(Matching_nomenclature_coding_pquery = str_squish(merged_df$c.nomenclature) == str_squish(merged_df$hgvsc_protein))

merged_df <- merged_df %>%
  mutate(Matching_nomenclature_protein_pquery  = str_squish(merged_df$p.nomenclature) == str_squish(merged_df$hgvsp_protein))

## Save this and manually inspect variants which do not match using Alamut. 
# Add a comments column. 
merged_df$`Comments on nomenclature` <- ""

# Add column to determine if to proceed. 
merged_df$`To proceed?` <- ""

# Add column to dictate what column to proceed with. 
merged_df$`Lab, coding, or protein?` <- ""

# Add comment to highlight splicing variants and so those which may have no available p.
merged_df$`Comments on nomenclature`[grepl("\\+|\\-", merged_df$c.nomenclature)] <- "Splicing defect."

# Save this as a csv to inspect. 
#write.csv(merged_df, "Results_variant_recoder_merged_inspect_df.csv", row.names = FALSE, quote = TRUE)

## Reload manually inspecting merged data frame. 
Results_variant_recoder_protein_df_EDITS <- read.csv("Results_variant_recoder_merged_inspect_df_EDIT.csv", nrows = 626)

## For variants where all the variant recoder (coding and protein) queries are the same as the lab
## reported variant, or where the protein query hasn't retrieved (NA) but coding variant recoder query 
## retrieved c. and p. matches the lab, proceed. 
with(Results_variant_recoder_protein_df_EDITS, {
     condition <- Matching_nomenclature_coding_cquery & 
     Matching_nomenclature_protein_cquery & 
    (is.na(Matching_nomenclature_coding_pquery) | Matching_nomenclature_coding_pquery) &
    (is.na(Matching_nomenclature_protein_pquery) | Matching_nomenclature_protein_pquery)
    
    # Only update rows where condition is TRUE and current value is NA
    idx <- condition & is.na(To.proceed.)
    Results_variant_recoder_protein_df_EDITS$To.proceed.[idx] <<- TRUE
    
    idx2 <- condition & is.na(Lab..coding..or.protein.)
    Results_variant_recoder_protein_df_EDITS$Lab..coding..or.protein.[idx2] <<- TRUE
    
    idx3 <- condition & Lab..coding..or.protein. == ""
    Results_variant_recoder_protein_df_EDITS$Lab..coding..or.protein.[idx3] <<- "Coding"
  })

## Check variants which failed checks. 
Failed_check <- Results_variant_recoder_protein_df_EDITS[Results_variant_recoder_protein_df_EDITS$To.proceed. == FALSE,]

## Check variants labelled "Other" in `Lab..coding..or.protein.` column as these require further checks.
Other <- Results_variant_recoder_protein_df_EDITS %>%
         filter(Lab..coding..or.protein. == "Other")

# Edit entry from Oxford with corrected nomenclature. 
Results_variant_recoder_protein_df_EDITS[Results_variant_recoder_protein_df_EDITS$ID == 440,
                                         c("Lab..coding..or.protein.", "hgvsg","hgvsg_alt", "hgvsc", "hgvsp")] <- list("Coding", "NC_000007.14:g.150957445_150957446insTAATGGTAGCGCA",
                                                                                                                       "NC_000007.14:g.150957436_150957437insGGTAGCGCATAAT", 
                                                                                                                       "NM_000238.4:c.982_983insATTATGCGCTACC", 
                                                                                                                       "NP_000229.1:p.Arg328HisfsTer8")

## Assign final HGVS to use for each variant. 
# Select only variants which were to proceed. 
Results_variant_recoder_protein_df_final <- Results_variant_recoder_protein_df_EDITS[Results_variant_recoder_protein_df_EDITS$To.proceed. == TRUE,]

# Create columns to hold final HGVS info.
Results_variant_recoder_protein_df_final$Final_gnomen <- NA
Results_variant_recoder_protein_df_final$Final_altgnomen <- NA
Results_variant_recoder_protein_df_final$Final_cnomen <- NA
Results_variant_recoder_protein_df_final$Final_pnomen <- NA

# For each source of data (lab, c nomen variant recoder query, p nomen variant recoder query) assign to final hgvs according to previous inspection of data.
Results_variant_recoder_protein_df_final[Results_variant_recoder_protein_df_final$Lab..coding..or.protein. == "Lab", c("Final_gnomen", "Final_cnomen","Final_pnomen")] <- Results_variant_recoder_protein_df_final[Results_variant_recoder_protein_df_final$Lab..coding..or.protein. == "Lab", 
                                                                                                                                                                                                                   c("g.nomenclature", "c.nomenclature", "p.nomenclature")]

Results_variant_recoder_protein_df_final[Results_variant_recoder_protein_df_final$Lab..coding..or.protein. == "Coding", c("Final_gnomen", "Final_altgnomen", "Final_cnomen","Final_pnomen")] <- Results_variant_recoder_protein_df_final[Results_variant_recoder_protein_df_final$Lab..coding..or.protein. == "Coding", 
                                                                                                                                                                                                                   c("hgvsg", "hgvsg_alt", "hgvsc", "hgvsp")]

Results_variant_recoder_protein_df_final[Results_variant_recoder_protein_df_final$Lab..coding..or.protein. == "Protein", c("Final_gnomen", "Final_altgnomen", "Final_cnomen","Final_pnomen")] <- Results_variant_recoder_protein_df_final[Results_variant_recoder_protein_df_final$Lab..coding..or.protein. == "Protein", 
                                                                                                                                                                                                                   c("hgvsg_protein", "hgvsg_alt_protein", "hgvsc_protein", "hgvsp_protein")]

## Save the final version of the data.
#write.csv(Results_variant_recoder_protein_df_final, "Results_variant_recoder_protein_df_final.csv", row.names = FALSE, quote = TRUE)

## Create a unique c. list for a later script which retrieves variant info from Alamut. 
Results_variant_recoder_protein_df_final <- read.csv("Results_variant_recoder_protein_df_final.csv")

c_list <- unique(Results_variant_recoder_protein_df_final$`Final_cnomen`)

c_list <- gsub("NM_000238\\.4:", "", c_list)

writeLines(as.character(c_list), "All_variants_cnomen.csv")

## Check the number of unique variants after HGVS filtering. 
############################################################

## Aberdeen.
Aberdeen_df_filtered <- Results_variant_recoder_protein_df_final[
  Results_variant_recoder_protein_df_final$Lab == "Aberdeen", 
]

# Determine number of unique variants. 
length(unique(Aberdeen_df_filtered$Final_cnomen))

Aberdeen_missense_variants_filtered <- Aberdeen_df_filtered[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Aberdeen_df_filtered$Final_pnomen
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Aberdeen_df_filtered$Final_pnomen
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Aberdeen_df_filtered$Final_cnomen
    ),
]

length(unique(Aberdeen_missense_variants_filtered$Final_cnomen))


## Belfast.
Belfast_df_filtered <- Results_variant_recoder_protein_df_final[
  Results_variant_recoder_protein_df_final$Lab == "Belfast", 
]

# Determine number of unique variants. 
length(unique(Belfast_df_filtered$Final_cnomen))

Belfast_missense_variants_filtered <- Belfast_df_filtered[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Belfast_df_filtered$Final_pnomen
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Belfast_df_filtered$Final_pnomen
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Belfast_df_filtered$Final_cnomen
    ),
]

length(unique(Belfast_missense_variants_filtered$Final_cnomen))


## Brompton.
Brompton_df_filtered <- Results_variant_recoder_protein_df_final[
  Results_variant_recoder_protein_df_final$Lab == "Brompton", 
]

# Determine number of unique variants. 
length(unique(Brompton_df_filtered$Final_cnomen))

Brompton_missense_variants_filtered <- Brompton_df_filtered[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Brompton_df_filtered$Final_pnomen
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Brompton_df_filtered$Final_pnomen
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Brompton_df_filtered$Final_cnomen
    ),
]

length(unique(Brompton_missense_variants_filtered$Final_cnomen))


## Manchester.
Manchester_df_filtered <- Results_variant_recoder_protein_df_final[
  Results_variant_recoder_protein_df_final$Lab == "Manchester", 
]

# Determine number of unique variants. 
length(unique(Manchester_df_filtered$Final_cnomen))

Manchester_missense_variants_filtered <- Manchester_df_filtered[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Manchester_df_filtered$Final_pnomen
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Manchester_df_filtered$Final_pnomen
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Manchester_df_filtered$Final_cnomen
    ),
]

length(unique(Manchester_missense_variants_filtered$Final_cnomen))


## Oxford.
Oxford_df_filtered <- Results_variant_recoder_protein_df_final[
  Results_variant_recoder_protein_df_final$Lab == "Oxford", 
]

# Determine number of unique variants. 
length(unique(Oxford_df_filtered$Final_cnomen))

Oxford_missense_variants_filtered <- Oxford_df_filtered[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Oxford_df_filtered$Final_pnomen
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Oxford_df_filtered$Final_pnomen
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Oxford_df_filtered$Final_cnomen
    ),
]

length(unique(Oxford_missense_variants_filtered$Final_cnomen))


## Total number of unique missense variants. 
length(unique(Results_variant_recoder_protein_df_final$Final_cnomen))

Total_missense_variants_filtered <- Results_variant_recoder_protein_df_final[
  grepl(
    "^NP_.*:p\\.([A-Z][a-z]{2}|[A-Z])[0-9]+([A-Z][a-z]{2}|[A-Z]|Ter|\\*)$",
    Results_variant_recoder_protein_df_final$Final_pnomen
  ) &
    !grepl(
      "^NP_.*:p\\.((([A-Z][a-z]{2})([0-9]+)\\3)|(([A-Z])([0-9]+)\\6)|(([A-Z][a-z]{2}|[A-Z])[0-9]+=)|.*Ter$|.*\\*$)",
      Results_variant_recoder_protein_df_final$Final_pnomen
    ) &
    grepl(
      "^NM_.*:c\\.[0-9]+[ACGT]>[ACGT]$",
      Results_variant_recoder_protein_df_final$Final_cnomen
    ),
]

length(unique(Total_missense_variants_filtered$Final_cnomen))

### Visualise the classifications for each variant
variant_classifications <- Total_missense_variants_filtered[, c(
  "Final_cnomen",
  "Lab.variant.classification",
  "Classification.date"
)]

# Remove rows that are completely identical across all three columns.
variant_classifications <- unique(variant_classifications)

# Merge classifications and dates for each unique variant.
variant_classifications_combined <- aggregate(
  cbind(Lab.variant.classification, Classification.date) ~ Final_cnomen,
  data = variant_classifications,
  FUN = function(x) paste(unique(x), collapse = "; ")
)

# Count the number of VUS rows.
sum(grepl("VUS", variant_classifications_combined$Lab.variant.classification))

## Create anonymised collated dataset file for supplementary figures.
# For the repository, this has been edited to not show the labs associated anon number. This will not run without editing.
Results_variant_recoder_protein_df_final$Lab <- dplyr::recode(
  Results_variant_recoder_protein_df_final$Lab,
  "A" = "Lab 1",
  "B" = "Lab 2",
  "C" = "Lab 3",
  "D" = "Lab 4",
  "E" = "Lab 5"
)

# Remove identifiable / unnecessary columns.
Results_variant_recoder_protein_df_final_anonymised <- 
  Results_variant_recoder_protein_df_final[, -c(2, 4, 10, 12, 13, 17:ncol(Results_variant_recoder_protein_df_final))]

# Save anonymised dataset.
write.csv(Results_variant_recoder_protein_df_final_anonymised, "Results_variant_recoder_protein_df_final_anonymised.csv", row.names = FALSE,quote = TRUE)
