% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregateMSexperiment.R
\name{aggregateProteins}
\alias{aggregateProteins}
\title{aggregate peptides - given msexperiment with transitons}
\usage{
aggregateProteins(msexp, FUN = sum)
}
\arguments{
\item{msexp}{- msexperiment}

\item{FUN}{aggregation function to use}
}
\description{
aggregate peptides - given msexperiment with transitons
}
\examples{
data( feature_alignment_requant )
msexp = loadTransitonsMSExperiment( feature_alignment_requant , nrt = 3 , peptop = 3 )
x = msexp
x = removeDecoys( x )
table( table( x$pepinfo$ProteinName ) )
y = aggregateProteins( x )
stopifnot(rownames(y$pepinfo$ProteinName) == rownames(y$Intensity))
stopifnot( table( table( y$pepinfo$ProteinName ) ) == 50 )

}
\seealso{
aggregatePeptide
}

