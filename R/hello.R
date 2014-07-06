#'A package with some utility functions for the analysis of swath data
#'@name imsbInfer
#'@docType package
#'@import roxygen2 data.table gplots reshape2
NULL
#' @name feature_alignment_requant
#' @title feature_alignment_requant - output of feature aligner requant
#' @description this is an example output when running feature aligner requant on swath data
#' @docType data
#' @usage feature_alignment_requant
#' @format see example below
#' @author Lorenz Blum
#' @seealso \code{\link{convertLF2Wideformat}} and \code{\link{read2msExperiment.data.frame}}
#' data(feature_alignment_requant)
#' colnames(feature_alignment_requant)
#' head(feature_alignment_requant[,1:8])
NULL
#' @name SDat
#' @title SDat - example of an msexperiment class
#' @description this is an example output of calling the function read2msExperiment
#' @docType data
#' @usage feature_alignment_requant
#' @format see example below
#' @author Witold Wolski
#' @seealso \code{\link{feature_alignment_requant}} and \code{\link{read2msExperiment.data.frame}}
#' @examples
#' data(feature_alignment_requant)
#' SDat = read2msExperiment(feature_alignment_requant)
#' names(SDat)
NULL
# hack to supress _no visible binding for global variable _ warning in R CMD check.
utils::globalVariables(c("transition_group_id","ProteinName","align_origfilename","Intensity","aggr_Fragment_Annotation"), add = TRUE)
Sys.setlocale("LC_COLLATE", "C")