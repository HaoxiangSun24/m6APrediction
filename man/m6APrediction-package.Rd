\name{m6APrediction-package}
\alias{m6APrediction-package}
\alias{m6APrediction}
\docType{package}
\title{m6APrediction: Predict N6-Methyladenosine RNA Modification Sites}
\description{
A random forest framework for predicting N6-methyladenosine (m6A) RNA
modification sites from transcript and sequence features.
}
\details{
The package provides two exported functions:
\itemize{
  \item \code{\link{prediction_single}} predicts a single candidate site.
  \item \code{\link{prediction_multiple}} predicts one or more sites supplied
  as a data frame.
}

The bundled model was trained on 48,055 human m6A sites using six transcript
or sequence features and the five nucleotides of the DNA 5-mer centred on the
candidate adenosine.
}
\author{
Your Name
}
\keyword{package}
