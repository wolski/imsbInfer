.OpenMSPrecursorDefsMapping <- list("Filename"="align_origfilename",
                    "ProteinName"="ProteinName",
                    "Decoy"="decoy",
                    "StrippedSequence"="pep_sequence",
                    "ModifiedSequence"="FullPeptideName",
                    #"IsotopeLabelType"=NULL,
                    "PrecursorCharge" = "pep_charge",
                    "PrecursorMZ"= "mz",
                    "PrecursorRT"="RT",
                    "PrecursorScore"="m_score")

.OpenMSFragmentDefsMapping <-c("FragmentIonType"="ion_type",
                  "FragmentCharge"="frag_charge",
                  "FragmentIntensity"="aggr_Peak_Area")

.prepareDF <- function (df)
{
  colnames(df) <- gsub("m/z","mz",colnames(df))
  colnames(df) <- tolower(colnames(df))
  required = tolower(c("ProteinName","transition_group_id", "align_origfilename","decoy","FullPeptideName",
                                  "RT", "mz", "Intensity", "ProteinName", "m_score", "aggr_Fragment_Annotation",
                                  "aggr_Peak_Area" ))
  x = match(required, colnames(df))
  if(sum(is.na(x)) > 0){
    warning("missing required columns : ", required[is.na(x)])
  }
  df = df[, x]
  df$align_origfilename <- gsub("_with_dscore_filtered.csv","", basename(df$align_origfilename))
  return(df)
}

#' Prepare data from OpenSwathOutput
#' @export
#' @examples 
#' rm(list=ls())
#' library(readr)
#' library(imsbInfer2)
#' data2 <- read_tsv("inst/extdata/example.tsv.gz",col_names = TRUE)
#' head(data2)
#' 
#' far <- .prepareDF(data2)
#' head(far)
#' tmp[1:3]
#' data3 <- prepareOpenSwathData(data2)
#' 
prepareOpenSwathData <- function(far){
  far <- .prepareDF(far)
  apa = as.character(far$aggr_peak_area)
  afa = as.character(far$aggr_fragment_annotation)
  # split transition intensities
  
  #transints = lapply(apa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  transints = strsplit(apa,";",fixed=TRUE)
  
  # split transition names
  #transids = lapply(afa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  transids = strsplit(afa,";",fixed=TRUE)  
  
  # Fix protein names
  protnames <-strsplit(far$proteinname, split="/", fixed=TRUE)
  protnames <-lapply(protnames, function(x){ c(x[1], sort(x[2:length(x)])) })
  protnames <- sapply(protnames, function(x){paste(x,collapse="/")})
  far$proteinname <- protnames
  
  # prepare output
  lx = length(transids)
  stopifnot(lx == nrow(far))

  far <- far[, !(names(far) %in% c("aggr_peak_area","aggr_fragment_annotation"))]

  # extend
  lengths <-sapply(transids,length)
  idx <- rep(1:lx, lengths)
  far <- far[idx, ]

  #
  transids <- unlist(transids)
  transids <- gsub("DECOY_","DECOY", transids, fixed=TRUE )
  cnamessplit <- strsplit(as.character(transids),split="_",fixed=TRUE)
  transids <- do.call("rbind",cnamessplit)
  colnames(transids)<-c("frag_id", "ion_type","frag_charge","pep_sequence","pep_charge")
  length(as.numeric(unlist(transints)))
  far <- data.frame(far, aggr_peak_area = as.numeric(unlist(transints)), transids )
  
  far <-far[,tolower(unlist(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping)))]
  names(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping))
  colnames(far) <- names(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping))
  far <- data.frame(far, IsotopeLabelType = "L")
  return(far)
}

