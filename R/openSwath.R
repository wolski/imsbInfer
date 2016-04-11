split2long <- function(far){
  far <- data2
  #far = feature_alignment_requant
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # prepare output
  lx = length(transids)
  stopifnot(lx == nrow(far))

  far<- far[, !(names(far) %in% c("aggr_Peak_Area","aggr_Fragment_Annotation"))]
  message("prepared dataframe")
  # extend
  lengths <-sapply(transids,length)
  idx <- rep(1:lx, lengths)
  range(idx)
  far <- far[idx, ]
  message("and finally nrow " , nrow(far), " ncol ", ncol(far))
  far <- data.frame(far, aggr_Peak_Area = as.numeric(unlist(transints)),
                    aggr_Fragment_Annotation =  as.character(unlist( transids) ) )
  return(far)
}

prepareDF <- function (df)
{
  colnames(df) <- gsub("m/z","mz",colnames(df))
  required = c("transition_group_id", "align_origfilename","decoy",
               "RT", "mz", "Intensity", "ProteinName", "m_score", "aggr_Fragment_Annotation",
               "aggr_Peak_Area" )
  x = match(tolower(required), tolower(colnames(df)))
  stopifnot(required == colnames(df)[x])
  df = df[, x]
  df <- df[order(df$transition_group_id), ]
  return(df)
}
