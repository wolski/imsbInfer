% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/swathDataUtils.R
\name{keepRTRange.msexperiment}
\alias{keepRTRange.msexperiment}
\title{filter data given RT range}
\usage{
\method{keepRTRange}{msexperiment}(obj, rtrange = c(1000, 7000), ...)
}
\arguments{
\item{rtrange}{rt range to keep default 1000 - 7000}

\item{data}{msexperiment}
}
\description{
filter data given RT range
}
\examples{
data(feature_alignment_requant)
SDat = read2msExperiment(feature_alignment_requant)
dim(SDat)
SD2 = keepRTRange(SDat,c(3000,4000))
dim(SD2)
}
