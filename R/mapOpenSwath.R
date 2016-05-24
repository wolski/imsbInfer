.OpenMSPrecursorDefsMapping <- list("FileName"="align_origfilename",
                    "ProteinName"="ProteinName",
                    "Decoy"="decoy",
                    "StrippedSequence"="Sequence",
                    "ModifiedSequence"="FullPeptideName",
                    "PrecursorCharge" = "pep_charge",
                    "PrecursorMZ"= "mz",
                    "PrecursorRT"="RT",
                    "PrecursorScore"="m_score",
                    "MS2IntensityAggregated" = "Intensity")

.OpenMSFragmentDefsMapping <-c("FragmentIonType"="ion_type",
                  "FragmentCharge"="frag_charge",
                  "FragmentIntensity"="aggr_Peak_Area")

.prepareDF <- function (df)
{
  colnames(df) <- gsub("m/z","mz",colnames(df))
  colnames(df) <- tolower(colnames(df))
  required = tolower(c("ProteinName","transition_group_id", "align_origfilename","decoy","FullPeptideName","Sequence",
                                  "RT", "mz", "Intensity", "m_score", "aggr_Fragment_Annotation",
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
#' data2 <- read_tsv(file.path(path.package("imsbInfer2"),"extdata/example.tsv.gz"),col_names = TRUE)
#' head(data2)
#' 
#' #far <- .prepareDF(data2)
#' #far %>% glimpse
#' #head(far)
#' #tmp[1:3]
#' data3 <- prepareOpenSwathData(data2)
#' colnames(data3)
prepareOpenSwathData <- function(far){
  far <- .prepareDF(far)

  message("selected essential columns")
  apa = as.character(far$aggr_peak_area)
  afa = as.character(far$aggr_fragment_annotation)
  # split transition intensities
  transints = strsplit(apa,";",fixed=TRUE)
  message("prepared transition intensity")
  # split transition names
  transids = strsplit(afa,";",fixed=TRUE)
  message("prepared transition ids")
  
  # Fix protein names
  protnames <-strsplit(far$proteinname, split="/", fixed=TRUE)
  protnames <-lapply(protnames, function(x){ c(x[1], sort(x[2:length(x)])) })
  protnames <- sapply(protnames, function(x){paste(x,collapse="/")})
  far$proteinname <- protnames
  message("adjusting protein names done")

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
  
  far <- data.frame(far, aggr_peak_area = as.numeric(unlist(transints)), transids )
  
  
  print(setdiff(tolower(unlist(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping))) , colnames(far)))
  print(setdiff(colnames(far),tolower(unlist(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping)))))
  
  far <-far[,tolower(unlist(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping)))]
  
  
  colnames(far) <- names(c(.OpenMSPrecursorDefsMapping,.OpenMSFragmentDefsMapping))
  far <- data.frame(far, LabelType = "L", FragmentInterference = NA, MS1Intensity = NA)
  return(far)
}

