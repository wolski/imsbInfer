.formalOpenSwathMap <-list("z" = "pep_charge",
                          "RT"="RT",
                          "ModifiedSequence"="FullPeptideName",
                          "Sequence"="pep_sequence",
                          "Filename"="align_origfilename",
                          "Decoy"="decoy", 
                          "mz"= "mz",
                          "Score"="m_score",
                          "IonType"="ion_type",
                          "FragZ"="frag_charge",
                          "Intensity"="aggr_Peak_Area")


prepareDF <- function (df)
{
  
  colnames(df) <- gsub("m/z","mz",colnames(df))
  required = c("transition_group_id", "align_origfilename","decoy","FullPeptideName",
               "RT", "mz", "Intensity", "ProteinName", "m_score", "aggr_Fragment_Annotation",
               "aggr_Peak_Area" )
  
  x = match(tolower(required), tolower(colnames(df)))
  stopifnot(required == colnames(df)[x])
  df = df[, x]
  df$align_origfilename <- gsub("_with_dscore_filtered.csv","", basename(df$align_origfilename))
  return(df)
}

prepareOpenSwathData <- function(far){
  far <- data2
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # prepare output
  lx = length(transids)
  stopifnot(lx == nrow(far))

  far <- far[, !(names(far) %in% c("aggr_Peak_Area","aggr_Fragment_Annotation"))]
  
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
  far <- data.frame(far, aggr_Peak_Area = as.numeric(unlist(transints)), transids )
  
  far <-far[,unlist(.formalOpenSwathMap)]
  colnames(far) <- names(.formalOpenSwathMap)
  
  return(far)
}

