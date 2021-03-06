% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swathDataUtils.R
\name{read2msExperiment.default}
\alias{read2msExperiment}
\alias{read2msExperiment.default}
\title{read 2 feature alginer long format and generate an msexperiment (takes peptide quantifications of feature aligner)}
\usage{
\method{read2msExperiment}{default}(obj, ...)
}
\arguments{
\item{obj}{filename to an openswath feature alignment output}

\item{...}{unused}
}
\value{
msexperiment
}
\description{
read 2 feature alginer long format and generate an msexperiment (takes peptide quantifications of feature aligner)
}
\examples{
data(feature_alignment_requant) 
SDat = read2msExperiment(feature_alignment_requant)
stopifnot(rownames(SDat$pepinfo)==rownames(SDat$Intensity))
\dontrun{res = read2msExperiment("path/to/feature_alignment_requant.tsv")}
}

