\name{dna_encoding}
\alias{dna_encoding}
\title{Encode DNA 5-mers as Nucleotide Factors}
\description{
Splits each five-nucleotide sequence into five factor columns
(\code{nt_pos1} to \code{nt_pos5}) with levels \code{A}, \code{T}, \code{C}
and \code{G}. These columns are the sequence features used by the bundled
random forest model.
}
\usage{
dna_encoding(dna_strings)
}
\arguments{
  \item{dna_strings}{Character vector of DNA 5-mers. Each element must contain
  five characters from A, T, C and G, and the third character must be
  \code{A} (the candidate adenosine). Uracil (U) is accepted and converted to
  T.}
}
\value{
A data frame with one row per input sequence and five factor columns.
}
\keyword{internal}
