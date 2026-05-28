### Set working directory.
##########################
setwd("~/KCNH2_project")

### Load libraries. 
###################
library(data.table)
library(stringr)
suppressWarnings(library('bedr', lib.loc = "/data/mukund_temp/AlphaMissense/"))

### Load in the KCNH2 lab variants.
###################################
## KCNH2 lab variant classifications. 
Results_variant_recoder_protein_df_final <- read.csv("Results_variant_recoder_protein_df_final.csv")
Results_variant_recoder_protein_df_final$Final_cnomen <- gsub("NM_000238\\.4:", "NM_000238:", Results_variant_recoder_protein_df_final$Final_cnomen)

Unique_variant_df <- Results_variant_recoder_protein_df_final[!duplicated(Results_variant_recoder_protein_df_final$Final_cnomen), 34:37]

### Set file paths for AlphaMissense predictions.
#################################################
am_dir       <- "/data/mukund_temp/AlphaMissense/"
am_pred_file <- file.path(am_dir, "AlphaMissense_hg38.tsv.gz")
am_gene_file <- file.path(am_dir, "AlphaMissense_gene_hg38.tsv.gz")

### Gather the computational pathogencity predictions. 
######################################################

### KEY NOTE. This HAS to be run on R lab server to work due to the VEP plugins used. 

## Create empty DF to hold the REVEL, SpliceAI, and AlphaMissense scores. 
Variants_computational_scores <- data.frame(Variant_cnomen = character(), REVEL_score = character(), Delta_Score_AG = character(), Delta_score_AL = character(), Delta_score_DG = character(),
                                   Delta_score_DL = character(), Delta_position_AG = character(), Delta_position_AL = character(), Delta_position_DG = character(), Delta_position_DL = character(),
                                   Transcript = character(), HGVSc = character(), HGVSp = character(), Chrom = character(), Pos = character(), Ref = character(), Alt = character(), 
                                   Transcript_ID = character(), Protein_variant = character(), AM_pathogenicity = character(), AM_class = character(), VEP_Comment = character(), AM_Comment = character())

## For loop to gather the computational scores for each variant using VEP and put this in the above dataframe. 
for (i in seq_len(nrow(Unique_variant_df))) {
  Unique_variant_df_row <- Unique_variant_df[i, ] 
  
  print(paste("Processing variant",Unique_variant_df_row$Final_cnomen))
  vep_cmd <- paste0("echo '", Unique_variant_df_row$Final_cnomen, 
                    "' | /data/resources/ensembl-vep/vep -i STDIN --format hgvs --cache --dir_cache /data/vep/ --dir_plugins /data/vep/plugins108/ -a GRCh38 --fa /data/humanGenome/GRCh38_no_alt/GRCh38_no_alt_analysis_set.fa --refseq --hgvs --symbol --tab --warning_file /tmp/vep_warnings.txt --force_overwrite --plugin REVEL,/data/vep/plugins108/data_REVEL/new_tabbed_revel_hg38.tsv.gz,1 --plugin SpliceAI,snv=/data/vep/plugins108/data_SpliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/data/vep/plugins108/data_SpliceAI/spliceai_scores.raw.indel.hg38.vcf.gz -o STDOUT")
  raw <- system(vep_cmd, intern = TRUE)
  raw2 <- as.data.frame(raw)
  
  vep_lines <- raw[grepl("^#Uploaded_variation|^NM_", raw)]
  
  vep_lines[1] <- sub("^#", "", vep_lines[1])
  
  vep_df <- read.delim(
  text = paste(vep_lines, collapse = "\n"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE)
  
  if(nrow(vep_df) == 0) {
    Row <- vep_df[NA,]
    
    Row$Delta_Score_AG <- NA
    Row$Delta_Score_AL <- NA
    Row$Delta_Score_DG <- NA
    Row$Delta_Score_DL <- NA
    
    Row$Delta_position_AG <- NA
    Row$Delta_position_AL <- NA
    Row$Delta_position_DG <- NA
    Row$Delta_position_DL <- NA
    
    Row$VEP_Comment <- "VEP search unsuccessful."
    
  } else {
  
    Row <- vep_df[vep_df$Feature == 'NM_000238.3',]
    
    splice_split <- strsplit(Row$SpliceAI_pred, "\\|")[[1]]
    
    Row$Delta_Score_AG <- splice_split[2]
    Row$Delta_Score_AL <- splice_split[3]
    Row$Delta_Score_DG <- splice_split[4]
    Row$Delta_Score_DL <- splice_split[5]
    
    Row$Delta_position_AG <- splice_split[6]
    Row$Delta_position_AL <- splice_split[7]
    Row$Delta_position_DG <- splice_split[8]
    Row$Delta_position_DL <- splice_split[9]
    
    Row$VEP_Comment <- ""
  }
  
    if (grepl("dup|del|ins", as.character(Unique_variant_df_row$Final_cnomen))) {
      am_result_final <- data.frame(CHROM = NA, POS = NA, REF = NA, ALT = NA, 
                                    transcript_id = NA, protein_variant = NA, 
                                    am_pathogenicity = NA, am_class = NA, AM_Comment = "No AM data, dup/del/ins variant.")
      
    } else {
    variant <- Unique_variant_df_row$Final_gnomen
    
    # Extract position
    position <- str_extract(variant, "(?<=g\\.)\\d+")
    
    # Extract ref allele
    ref <- str_extract(variant, "(?<=\\d)[ACGT]+")
    
    # Extract alt allele
    alt <- str_extract(variant, "(?<=\\>)[ACGT]+")
    
    # Build chr:pos format
    variant_formatted <- paste0("7:", position)
    
    convertBed <- function(variant_formatted) {
      out <- c()
      for (v in variant_formatted) {
        chromosome <- strsplit(v, ":")[[1]][1]
        position   <- strsplit(v, ":")[[1]][2]
        start      <- as.numeric(strsplit(position, "-")[[1]][1])
        end        <- strsplit(position, "-")[[1]][2]
        end        <- ifelse(is.na(end), start, as.numeric(end))
        start      <- start - 1
        if (!startsWith(chromosome, "chr")) 
          chromosome <- paste0("chr", chromosome)
          out <- c(out, paste0(chromosome, ":", start, "-", end))
      }
      return(out)
    }
    
    bed_query <- convertBed(variant_formatted)
    
    am_result <- tabix(bed_query, am_pred_file,
                       check.chr = FALSE, check.zero.based = FALSE,
                       deleteTmpDir = TRUE)
    
    if(is.null(am_result)){
      am_result_final <- data.frame(CHROM = NA, POS = NA, REF = NA, ALT = NA, 
                                    transcript_id = NA, protein_variant = NA, 
                                    am_pathogenicity = NA, am_class = NA, AM_Comment = "No AM data found.")
      
    } else {
      am_result_final <- am_result[am_result$REF == ref & am_result$ALT == alt,]
    }
    
    
    
    if(nrow(am_result_final) == 0){
      am_result_final <- am_result_final[NA,]
      am_result_final$AM_Comment = "No matching ref/alt."
    } 
    
    if(isTRUE(am_result_final$REF == ref) && isTRUE(am_result_final$ALT == alt)){
      am_result_final$AM_Comment = ""
    }
    }

  Variants_computational_scores <- rbind(Variants_computational_scores, data.frame(Variant_cnomen = Unique_variant_df_row$Final_cnomen, REVEL_score = Row$REVEL, Delta_Score_AG = Row$Delta_Score_AG, Delta_score_AL = Row$Delta_Score_AL, Delta_score_DG = Row$Delta_Score_DG,
                                                                 Delta_score_DL = Row$Delta_Score_DL, Delta_position_AG = Row$Delta_position_AG, Delta_position_AL = Row$Delta_position_AL, Delta_position_DG = Row$Delta_position_DG, 
                                                                 Delta_position_DL = Row$Delta_position_DL, Transcript = Row$Feature, HGVSc = Row$HGVSc, HGVSp = Row$HGVSp, Chrom = am_result_final$CHROM, Pos = am_result_final$POS, 
                                                                 Ref = am_result_final$REF, Alt = am_result_final$ALT, Transcript_ID = am_result_final$transcript_id, Protein_variant = am_result_final$protein_variant, 
                                                                 AM_pathogenicity = am_result_final$am_pathogenicity, AM_class = am_result_final$am_class, VEP_Comment = Row$VEP_Comment, AM_Comment = am_result_final$AM_Comment))
}

## Save the result. 
write.csv(Variants_computational_scores, "Variants_computational_scores.csv", row.names = FALSE, quote = TRUE)

### Determine REVEL evidence strengths (PP3/BP4).
#################################################

## Reload in the variant computational scores. 
Variants_computational_scores <- read.csv("Variants_computational_scores.csv")

## Define a dataframe to store the classification results for REVEL scores.
REVEL_evidence_result <- data.frame(Variant_cnomen = character(), 
                                    HGVSc = character(), 
                                    HGVSp = character(),
                                    REVEL_score = numeric(),
                                    PP3 = logical(), 
                                    BP4 = logical())

## Loop through variants and determine criteria based on REVEL thresholds.
for (i in seq_len(nrow(Variants_computational_scores))) {
  
  Variant_row <- Variants_computational_scores[i, ]
  
  Variant_cnomen <- Variant_row$Variant_cnomen
  HGVSc          <- Variant_row$HGVSc
  HGVSp          <- Variant_row$HGVSp
  REVEL_score    <- Variant_row$REVEL_score
  
  # Initialize logical flags for criteria
  PP3 <- FALSE
  BP4 <- FALSE
  
  # Categorize evidence strength according to REVEL thresholds
  if (!is.na(REVEL_score)) {
    if (REVEL_score >= 0.70) {
      PP3 <- TRUE
    } else if (REVEL_score <= 0.40) {
      BP4 <- TRUE
    }
  }
  
  # Create a final result row and append to the cumulative dataframe
  final_row <- data.frame(Variant_cnomen = Variant_cnomen, 
                          HGVSc = HGVSc, 
                          HGVSp = HGVSp,
                          REVEL_score = REVEL_score,
                          PP3 = PP3, 
                          BP4 = BP4)
  
  REVEL_evidence_result <- rbind(REVEL_evidence_result, final_row)
}

## Save this result.
write.csv(REVEL_evidence_result, "REVEL_evidence_result_KCNH2.csv", row.names = FALSE, quote = TRUE)


### Determine SpliceAI evidence strengths (PP3/BP4).
###################################################

## Define a dataframe to store the classification results for SpliceAI scores.
SpliceAI_evidence_result <- data.frame(Variant_cnomen = character(), 
                                       HGVSc = character(), 
                                       HGVSp = character(),
                                       Delta_Score_AG = numeric(),
                                       Delta_score_AL = numeric(),
                                       Delta_score_DG = numeric(),
                                       Delta_score_DL = numeric(),
                                       Delta_position_AG = numeric(),
                                       Delta_position_AL = numeric(),
                                       Delta_position_DG = numeric(),
                                       Delta_position_DL = numeric(),
                                       Max_Delta_Score = numeric(),
                                       PP3 = logical(), 
                                       BP4 = logical())

## Loop through variants and determine criteria based on SpliceAI thresholds.
for (i in seq_len(nrow(Variants_computational_scores))) {
  
  Variant_row <- Variants_computational_scores[i, ]
  
  Variant_cnomen    <- Variant_row$Variant_cnomen
  HGVSc             <- Variant_row$HGVSc
  HGVSp             <- Variant_row$HGVSp
  
  # Extract individual SpliceAI delta scores and positions
  Delta_Score_AG    <- Variant_row$Delta_Score_AG
  Delta_score_AL    <- Variant_row$Delta_score_AL
  Delta_score_DG    <- Variant_row$Delta_score_DG
  Delta_score_DL    <- Variant_row$Delta_score_DL
  
  Delta_position_AG <- Variant_row$Delta_position_AG
  Delta_position_AL <- Variant_row$Delta_position_AL
  Delta_position_DG <- Variant_row$Delta_position_DG
  Delta_position_DL <- Variant_row$Delta_position_DL
  
  # Initialize logical flags for criteria
  PP3 <- FALSE
  BP4 <- FALSE
  
  # Check if all 4 delta scores are missing to prevent the max() -Inf warning
  if (is.na(Delta_Score_AG) & is.na(Delta_score_AL) & is.na(Delta_score_DG) & is.na(Delta_score_DL)) {
    Max_Delta_Score <- NA
  } else {
    # Calculate the maximum delta score across all four splicing metrics safely
    Max_Delta_Score <- max(c(Delta_Score_AG, Delta_score_AL, Delta_score_DG, Delta_score_DL), na.rm = TRUE)
  }
  
  # Categorize evidence strength according to SpliceAI thresholds
  if (!is.na(Max_Delta_Score)) {
    if (Max_Delta_Score > 0.20) {
      PP3 <- TRUE
    } else if (Max_Delta_Score < 0.10) {
      BP4 <- TRUE
    }
  }
  
  # Create a final result row and append to the cumulative dataframe
  final_row <- data.frame(Variant_cnomen = Variant_cnomen, 
                          HGVSc = HGVSc, 
                          HGVSp = HGVSp,
                          Delta_Score_AG = Delta_Score_AG,
                          Delta_score_AL = Delta_score_AL,
                          Delta_score_DG = Delta_score_DG,
                          Delta_score_DL = Delta_score_DL,
                          Delta_position_AG = Delta_position_AG,
                          Delta_position_AL = Delta_position_AL,
                          Delta_position_DG = Delta_position_DG,
                          Delta_position_DL = Delta_position_DL,
                          Max_Delta_Score = Max_Delta_Score,
                          PP3 = PP3, 
                          BP4 = BP4)
  
  SpliceAI_evidence_result <- rbind(SpliceAI_evidence_result, final_row)
}

## Save this result.
write.csv(SpliceAI_evidence_result, "SpliceAI_evidence_result_KCNH2.csv", row.names = FALSE, quote = TRUE)