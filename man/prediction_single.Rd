\name{prediction_single}
\alias{prediction_single}
\title{Predict m6A Modification for a Single Site}
\description{
Predicts whether one candidate adenosine is an N6-methyladenosine (m6A) site
using a fitted random forest model. This is a convenient wrapper around
\code{\link{prediction_multiple}} for a single site.
}
\usage{
prediction_single(ml_fit = NULL, gc_content, RNA_type, RNA_region,
  exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer,
  positive_threshold = 0.5)
}
\arguments{
  \item{ml_fit}{A fitted \code{randomForest} object. If \code{NULL}
  (default), the model bundled with the package is used.}

  \item{gc_content}{Numeric scalar. GC content of the local sequence, between
  0 and 1.}

  \item{RNA_type}{Character scalar. Transcript biotype. One of
  \code{"mRNA"}, \code{"lincRNA"}, \code{"lncRNA"} or \code{"pseudogene"}.}

  \item{RNA_region}{Character scalar. Transcript region. One of
  \code{"CDS"}, \code{"intron"}, \code{"3'UTR"} or \code{"5'UTR"}.}

  \item{exon_length}{Numeric scalar. Length of the host exon, in the
  log-transformed units used during model training.}

  \item{distance_to_junction}{Numeric scalar. Distance from the site to the
  nearest splice junction, in the log-transformed units used during model
  training.}

  \item{evolutionary_conservation}{Numeric scalar. Cross-species conservation
  score, between 0 and 1.}

  \item{DNA_5mer}{Character scalar. Five nucleotides centred on the candidate
  adenosine, using A, T, C and G. The third nucleotide must be A. Uracil (U)
  is accepted and converted to T.}

  \item{positive_threshold}{Numeric scalar between 0 and 1. Probability cutoff
  used to assign the \code{"Positive"} class. Defaults to 0.5.}
}
\value{
A named character vector with elements \code{predicted_m6A_prob} and
\code{predicted_m6A_status}.
}
\details{
The bundled model is a random forest with 500 trees trained on 48,055 human
m6A sites. It uses six transcript or sequence features plus the five
nucleotide positions of the DNA 5-mer. Sites whose predicted probability
exceeds \code{positive_threshold} are reported as \code{"Positive"}.
}
\examples{
prediction_single(
  gc_content = 0.65,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 10.63,
  distance_to_junction = 7.66,
  evolutionary_conservation = 0.18,
  DNA_5mer = "GGACA"
)
}
\seealso{
\code{\link{prediction_multiple}} for predicting multiple sites.
}
\keyword{models}
\keyword{classif}
