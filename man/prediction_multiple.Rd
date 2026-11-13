\name{prediction_multiple}
\alias{prediction_multiple}
\title{Predict m6A Modification for Multiple Sites}
\description{
Applies a fitted random forest model to one or more candidate m6A sites. Each
site is described by six transcript or sequence features and a DNA 5-mer
centred on the candidate adenosine.
}
\usage{
prediction_multiple(ml_fit = NULL, feature_df, positive_threshold = 0.5)
}
\arguments{
  \item{ml_fit}{A fitted \code{randomForest} object. If \code{NULL}
  (default), the model bundled with the package is used.}

  \item{feature_df}{A data frame with the required feature columns.}

  \item{positive_threshold}{Numeric scalar between 0 and 1. Probability cutoff
  used to assign the \code{"Positive"} class. Defaults to 0.5.}
}
\value{
\code{feature_df} with two added columns:
\describe{
  \item{\code{predicted_m6A_prob}}{Predicted probability that each site is
  m6A-modified.}
  \item{\code{predicted_m6A_status}}{\code{"Positive"} or \code{"Negative"},
  assigned using \code{positive_threshold}.}
}
}
\details{
The required columns are \code{gc_content}, \code{RNA_type},
\code{RNA_region}, \code{exon_length}, \code{distance_to_junction},
\code{evolutionary_conservation} and \code{DNA_5mer}.

\code{RNA_type} must be one of \code{"mRNA"}, \code{"lincRNA"},
\code{"lncRNA"} or \code{"pseudogene"}. \code{RNA_region} must be one of
\code{"CDS"}, \code{"intron"}, \code{"3'UTR"} or \code{"5'UTR"}.
\code{DNA_5mer} must contain five nucleotides from A, T, C and G, with A at
the third position; uracil (U) is accepted and converted to T.
}
\examples{
sites <- data.frame(
  gc_content = c(0.65, 0.52),
  RNA_type = c("mRNA", "lncRNA"),
  RNA_region = c("CDS", "3'UTR"),
  exon_length = c(10.63, 9.13),
  distance_to_junction = c(7.66, 7.59),
  evolutionary_conservation = c(0.18, 0.62),
  DNA_5mer = c("GGACA", "GAACA"),
  stringsAsFactors = FALSE
)

prediction_multiple(feature_df = sites)
}
\seealso{
\code{\link{prediction_single}} for predicting a single site.
}
\keyword{models}
\keyword{classif}
