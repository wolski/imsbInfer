.OpenMSPrecursorDefsMapping <- list("Filename"="align_origfilename",
                    "ProteinName"="ProteinName",
                    "Decoy"="decoy",
                    "Sequence"="pep_sequence",
                    "ModifiedSequence"="FullPeptideName",
                    "IsotopeLabelType"=NULL,
                    "PrecursorCharge" = "pep_charge",
                    "PrecursorMZ"= "mz",
                    "PrecursorRT"="RT",
                    "PrecursorScore"="m_score")

.OpenMSFragmentDefsMapping <-c("FragmentIonType"="ion_type",
                  "FragmentCharge"="frag_charge",
                  "FragmentIntensity"="aggr_Peak_Area")

.prepareDF <- function (df)
{
  df<-data
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
#' library(readr)
#' library(imsbInfer2)
#' file = "d:/GoogleDrive/tissuecomparison/OpenSWATH/BAT_19strains/data/E1603291025_feature_alignment_requant.tsv.gz"
#' data <- read_tsv(file,col_names = TRUE)
#' data2 <- prepareOpenSwathData(data)
#' data <- .prepareDF(data)
#' 
prepareOpenSwathData <- function(far){
  data2 <- .prepareDF(far)
  far <- data2
  apa = as.character(far$aggr_peak_area)
  afa = as.character(far$aggr_fragment_annotation)
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # prepare output
  lx = length(transids)
  stopifnot(lx == nrow(far))
  
  far <- far[, !(names(far) %in% c("aggr_peak_Area","aggr_fragment_annotation"))]
  
  message("prepared dataframe")
  # extend
  lengths <-sapply(transids,length)
  idx <- rep(1:lx, lengths)
  range(idx)
  far <- far[idx, ]
  message("and finally nrow " , nrow(far), " ncol ", ncol(far))
  transids <- unlist(transids)
  transids <- gsub("DECOY_","DECOY", transids, fixed=TRUE )
  
  cnamessplit <- strsplit(as.character(transids),split="_",fixed=TRUE)
  transids <- do.call("rbind",cnamessplit)
  
  colnames(transids)<-c("frag_id", "ion_type","frag_charge","pep_sequence","pep_charge")
  far <- data.frame(far, aggr_peak_area = as.numeric(unlist(transints)), transids )
  
  far <-far[,unlist()]
  colnames(far) <- names()
  
  return(far)
}

