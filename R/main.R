

.m6a_cache <- new.env(parent = emptyenv())

.get_model <- function() {
  if (!is.null(.m6a_cache$fit)) {
    return(.m6a_cache$fit)
  }

  path <- system.file(
    "extdata", "rf_fit.rds",
    package = "m6APrediction", mustWork = TRUE
  )
  fit <- readRDS(path)
  .m6a_cache$fit <- fit
  fit
}


.stop <- function(...) {
  stop(..., call. = FALSE)
}

.check_threshold <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 || x > 1) {
    .stop("`positive_threshold` must be a single number between 0 and 1.")
  }
  as.numeric(x)
}

.check_numeric <- function(x, name, n) {
  if (!is.numeric(x)) {
    .stop(sprintf("`%s` must be numeric.", name))
  }
  if (length(x) != n || anyNA(x)) {
    .stop(sprintf("`%s` must contain %d non-missing value(s).", name, n))
  }
  as.numeric(x)
}

.check_factor <- function(x, name, levels) {
  if (is.factor(x)) {
    x <- as.character(x)
  }
  x <- trimws(as.character(x))

  if (anyNA(x)) {
    .stop(sprintf("`%s` contains missing values.", name))
  }

  bad <- setdiff(unique(x), levels)
  if (length(bad) > 0L) {
    .stop(sprintf(
      "Invalid value(s) in `%s`: %s. Allowed values: %s.",
      name,
      paste(bad, collapse = ", "),
      paste(levels, collapse = ", ")
    ))
  }

  factor(x, levels = levels)
}


#' Encode DNA 5-mers as nucleotide factors
#'
#' Splits each five-nucleotide sequence into five factor columns (`nt_pos1` to
#' `nt_pos5`) with levels `A`, `T`, `C` and `G`. These columns are the sequence
#' features used by the bundled random forest model.
#'
#' @param dna_strings Character vector of DNA 5-mers. Each element must contain
#'   five characters from A, T, C and G, and the third character must be `A`
#'   (the candidate adenosine). Uracil (U) is accepted and converted to T.
#'
#' @return A data frame with one row per input sequence and five factor columns.
#'
#' @keywords internal
dna_encoding <- function(dna_strings) {
  dna_strings <- toupper(trimws(as.character(dna_strings)))
  dna_strings <- gsub("U", "T", dna_strings, fixed = TRUE)

  if (length(dna_strings) == 0L) {
    .stop("`dna_strings` must contain at least one sequence.")
  }
  if (anyNA(dna_strings)) {
    .stop("`dna_strings` contains missing values.")
  }

  bad <- which(!grepl("^[ATCG]{5}$", dna_strings))
  if (length(bad) > 0L) {
    .stop(sprintf(
      paste0(
        "`dna_strings` must contain exactly five characters from A, T, C, G ",
        "(U is accepted and converted to T). Invalid row(s): %s."
      ),
      paste(utils::head(bad, 5L), collapse = ", ")
    ))
  }

  bad_center <- which(substr(dna_strings, 3L, 3L) != "A")
  if (length(bad_center) > 0L) {
    .stop(sprintf(
      paste0(
        "The third nucleotide of each 5-mer must be A ",
        "(the candidate adenosine). Invalid row(s): %s."
      ),
      paste(utils::head(bad_center, 5L), collapse = ", ")
    ))
  }

  nn <- nchar(dna_strings[1L])
  seq_m <- matrix(
    unlist(strsplit(dna_strings, "", fixed = TRUE)),
    ncol = nn, byrow = TRUE
  )
  colnames(seq_m) <- paste0("nt_pos", seq_len(nn))

  seq_df <- as.data.frame(seq_m, stringsAsFactors = FALSE)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  seq_df
}
-

.prepare_features <- function(feature_df) {
  if (!is.data.frame(feature_df)) {
    .stop("`feature_df` must be a data frame.")
  }
  if (nrow(feature_df) == 0L) {
    .stop("`feature_df` must contain at least one site.")
  }

  required <- c(
    "gc_content", "RNA_type", "RNA_region", "exon_length",
    "distance_to_junction", "evolutionary_conservation", "DNA_5mer"
  )
  missing_cols <- setdiff(required, names(feature_df))
  if (length(missing_cols) > 0L) {
    .stop(sprintf(
      "Missing required column(s): %s.",
      paste(missing_cols, collapse = ", ")
    ))
  }

  n <- nrow(feature_df)

  feature_df$gc_content <- .check_numeric(
    feature_df$gc_content, "gc_content", n
  )
  feature_df$exon_length <- .check_numeric(
    feature_df$exon_length, "exon_length", n
  )
  feature_df$distance_to_junction <- .check_numeric(
    feature_df$distance_to_junction, "distance_to_junction", n
  )
  feature_df$evolutionary_conservation <- .check_numeric(
    feature_df$evolutionary_conservation, "evolutionary_conservation", n
  )
  feature_df$RNA_type <- .check_factor(
    feature_df$RNA_type, "RNA_type",
    c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )
  feature_df$RNA_region <- .check_factor(
    feature_df$RNA_region, "RNA_region",
    c("CDS", "intron", "3'UTR", "5'UTR")
  )

  encoded <- dna_encoding(feature_df$DNA_5mer)

  model_df <- cbind(
    feature_df[, c(
      "gc_content", "RNA_type", "RNA_region", "exon_length",
      "distance_to_junction", "evolutionary_conservation"
    )],
    encoded
  )

  list(feature_df = feature_df, model_df = model_df)
}


#' Predict m6A modification for multiple sites
#'
#' Applies a fitted random forest model to one or more candidate m6A sites.
#' Each site is described by six transcript or sequence features and a DNA
#' 5-mer centred on the candidate adenosine.
#'
#' @param ml_fit A fitted `randomForest` object. If `NULL` (default), the model
#'   bundled with the package is used.
#' @param feature_df A data frame with the required feature columns.
#' @param positive_threshold Numeric scalar between 0 and 1. Probability cutoff
#'   used to assign the `"Positive"` class. Defaults to 0.5.
#'
#' @return `feature_df` with two added columns: `predicted_m6A_prob`
#'   (probability of being m6A-modified) and `predicted_m6A_status`
#'   (`"Positive"` or `"Negative"`).
#'
#' @details
#' The required columns are `gc_content`, `RNA_type`, `RNA_region`,
#' `exon_length`, `distance_to_junction`, `evolutionary_conservation` and
#' `DNA_5mer`.
#'
#' `RNA_type` must be one of `"mRNA"`, `"lincRNA"`, `"lncRNA"` or
#' `"pseudogene"`. `RNA_region` must be one of `"CDS"`, `"intron"`, `"3'UTR"`
#' or `"5'UTR"`. `DNA_5mer` must contain five nucleotides from A, T, C and G,
#' with A at the third position; uracil (U) is accepted and converted to T.
#'
#' @seealso [prediction_single()] for a single site.
#'
#' @examples
#' sites <- data.frame(
#'   gc_content = c(0.65, 0.52),
#'   RNA_type = c("mRNA", "lncRNA"),
#'   RNA_region = c("CDS", "3'UTR"),
#'   exon_length = c(10.63, 9.13),
#'   distance_to_junction = c(7.66, 7.59),
#'   evolutionary_conservation = c(0.18, 0.62),
#'   DNA_5mer = c("GGACA", "GAACA"),
#'   stringsAsFactors = FALSE
#' )
#'
#' prediction_multiple(feature_df = sites)
#'
#' @export
prediction_multiple <- function(ml_fit = NULL, feature_df,
                                positive_threshold = 0.5) {
  threshold <- .check_threshold(positive_threshold)

  prepared <- .prepare_features(feature_df)
  feature_df <- prepared$feature_df
  model_df <- prepared$model_df

  if (is.null(ml_fit)) {
    ml_fit <- .get_model()
  }
  if (!inherits(ml_fit, "randomForest")) {
    .stop("`ml_fit` must be a fitted randomForest object.")
  }

  prob <- predict(ml_fit, newdata = model_df, type = "prob")[, "Positive"]

  feature_df$predicted_m6A_prob <- as.numeric(prob)
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > threshold, "Positive", "Negative"
  )

  feature_df
}


#' Predict m6A modification for a single site
#'
#' Predicts whether one candidate adenosine is an N6-methyladenosine (m6A)
#' site using a fitted random forest model. This is a convenient wrapper around
#' [prediction_multiple()] for one site.
#'
#' @param ml_fit A fitted `randomForest` object. If `NULL` (default), the model
#'   bundled with the package is used.
#' @param gc_content Numeric scalar. GC content of the local sequence, between
#'   0 and 1.
#' @param RNA_type Character scalar. Transcript biotype. One of `"mRNA"`,
#'   `"lincRNA"`, `"lncRNA"` or `"pseudogene"`.
#' @param RNA_region Character scalar. Transcript region. One of `"CDS"`,
#'   `"intron"`, `"3'UTR"` or `"5'UTR"`.
#' @param exon_length Numeric scalar. Length of the host exon, in the
#'   log-transformed units used during model training.
#' @param distance_to_junction Numeric scalar. Distance from the site to the
#'   nearest splice junction, in the log-transformed units used during model
#'   training.
#' @param evolutionary_conservation Numeric scalar. Cross-species conservation
#'   score, between 0 and 1.
#' @param DNA_5mer Character scalar. Five nucleotides centred on the candidate
#'   adenosine, using A, T, C and G. The third nucleotide must be A. Uracil (U)
#'   is accepted and converted to T.
#' @param positive_threshold Numeric scalar between 0 and 1. Probability cutoff
#'   used to assign the `"Positive"` class. Defaults to 0.5.
#'
#' @return A named character vector with `predicted_m6A_prob` and
#'   `predicted_m6A_status`.
#'
#' @seealso [prediction_multiple()] for multiple sites.
#'
#' @examples
#' prediction_single(
#'   gc_content = 0.65,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 10.63,
#'   distance_to_junction = 7.66,
#'   evolutionary_conservation = 0.18,
#'   DNA_5mer = "GGACA"
#' )
#'
#' @export
prediction_single <- function(ml_fit = NULL, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  values <- list(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )

  lengths <- vapply(values, length, integer(1L))
  if (any(lengths != 1L)) {
    .stop(
      paste0(
        "prediction_single() accepts one site at a time. ",
        "Use prediction_multiple() to predict multiple sites."
      )
    )
  }

  feature_df <- as.data.frame(values, stringsAsFactors = FALSE)

  result <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = feature_df,
    positive_threshold = positive_threshold
  )

  c(
    predicted_m6A_prob = unname(result$predicted_m6A_prob),
    predicted_m6A_status = unname(result$predicted_m6A_status)
  )
}
