% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swathDataUtils.R
\name{setRequantToNA.msexperiment}
\alias{setRequantToNA.msexperiment}
\title{removes requant values from dataset}
\usage{
\method{setRequantToNA}{msexperiment}(obj)
}
\arguments{
\item{obj}{experiment}
}
\description{
removes requant values from dataset
}
\examples{
data(SDat)
SDatRT=orderByRT(SDat)
SDatRT=setRequantToNA(SDatRT)
SDatRT$RT[1:10] 
SDat = setRequantToNA(SDat)
head(SDat$rt)
}

