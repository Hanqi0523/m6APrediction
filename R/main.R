#' Encode DNA strings to positional nucleotide factors
#'
#' @description
#' Convert equal-length DNA strings (A/T/C/G) into a data.frame where each
#' column is a nucleotide position (`nt_pos1`, `nt_pos2`, …) and each column is
#' a factor with levels `c("A","T","C","G")`.
#'
#' @param dna_strings Character vector. All elements must have the **same**
#'   length and contain only A/T/C/G.
#'
#' @return A `data.frame` of factors with columns `nt_pos1 ... nt_posN` where
#'   `N = nchar(dna_strings[1])`. Each column has levels `c("A","T","C","G")`.
#'
#' @details
#' This helper is used by downstream prediction functions to transform the
#' 5-mer sequence into categorical predictors expected by the model.
#'
#' @examples
#' dna_encoding(c("ATCGA", "TTTTT"))
#' @export
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Batch prediction of m6A sites from a feature table
#'
#' @description
#' Use a fitted classifier and a feature table to compute m6A probability
#' and label per row. `DNA_5mer` will be expanded into `nt_pos1..nt_pos5`.
#'
#' @param ml_fit A fitted classifier supporting `predict(..., type = "prob")`.
#' @param feature_df A data.frame with columns:
#'   `gc_content`, `RNA_type`, `RNA_region`, `exon_length`,
#'   `distance_to_junction`, `evolutionary_conservation`, `DNA_5mer`.
#' @param positive_threshold Numeric in [0,1], default 0.5.
#' @return The input `feature_df` with two extra columns:
#'   `predicted_m6A_prob` and `predicted_m6A_status`.
#'
#' @details
#' `RNA_type` levels: `mRNA`, `lincRNA`, `lncRNA`, `pseudogene`;
#' `RNA_region` levels: `CDS`, `intron`, `3'UTR`, `5'UTR`.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' # ---- Load model & example CSV shipped in inst/extdata/ ----
#' rf_path  <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' csv_path <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")
#' ml_fit   <- readRDS(rf_path)
#' df       <- read.csv(csv_path, check.names = FALSE)
#'
#' # ---- Run batch prediction (copy/paste to console should work) ----
#' res <- prediction_multiple(ml_fit, feature_df = df, positive_threshold = 0.5)
#' head(res)
#' table(res$predicted_m6A_status)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame

  feature_df$RNA_type   <- factor(feature_df$RNA_type,
                                  levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  dna_df <- dna_encoding(feature_df$DNA_5mer)  # 会把 nt_pos1-5 设为因子，levels=c("A","T","C","G")
  for(i in 1:5){
    nm <- paste0("nt_pos", i)
    dna_df[[nm]] <- factor(dna_df[[nm]], levels = c("A","T","C","G"))
  }
  pred_input <- cbind(feature_df, dna_df)

  pred_cols <- c("gc_content", "RNA_type", "RNA_region",
                 "exon_length", "distance_to_junction",
                 "evolutionary_conservation",
                 "nt_pos1", "nt_pos2", "nt_pos4", "nt_pos5")

  prob_mat <- predict(ml_fit, newdata = pred_input[, pred_cols], type = "prob")
  pos_prob <- if ("Positive" %in% colnames(prob_mat)) prob_mat[, "Positive"] else prob_mat[, 2]
  pred_lab <- ifelse(pos_prob > positive_threshold, "Positive", "Negative")

  feature_df$predicted_m6A_prob   <- pos_prob
  feature_df$predicted_m6A_status <- factor(pred_lab, levels = c("Negative", "Positive"))

  return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}


#' Single-sample m6A prediction
#'
#' @description
#' A convenience wrapper that builds a one-row feature table and calls
#' [prediction_multiple()].
#'
#' @param ml_fit Fitted classifier as above.
#' @param gc_content Numeric.
#' @param RNA_type One of `mRNA`, `lincRNA`, `lncRNA`, `pseudogene`.
#' @param RNA_region One of `CDS`, `intron`, `3'UTR`, `5'UTR`.
#' @param exon_length Numeric.
#' @param distance_to_junction Numeric.
#' @param evolutionary_conservation Numeric.
#' @param DNA_5mer Character 5-mer (A/T/C/G).
#' @param positive_threshold Numeric in [0,1], default 0.5.
#' @return Named vector of length 2: `predicted_m6A_prob`, `predicted_m6A_status`.
#' @export
#'
#' @examples
#' # Load model shipped with the package
#' rf_path <- system.file("extdata", "rf_fit.rds", package = "m6APrediction")
#' ml_fit  <- readRDS(rf_path)
#'
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.55,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 1200,
#'   distance_to_junction = 45,
#'   evolutionary_conservation = 0.82,
#'   DNA_5mer = "ATCGA",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){

  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  pred_df <- prediction_multiple(ml_fit, feature_df, positive_threshold)

  returned_vector <- c(
    predicted_m6A_prob   = pred_df$predicted_m6A_prob[1],
    predicted_m6A_status = as.character(pred_df$predicted_m6A_status[1])
  )
  return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}


