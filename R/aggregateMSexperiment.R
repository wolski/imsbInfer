#' aggregate peptides - given msexperiment with transitons
#' @export
#' @param msexp msexperiment object
#' @param FUN aggregation function
#' @examples
#' data(feature_alignment_requant)
#' dim(feature_alignment_requant)/3
#' msexp = loadTransitonsMSExperiment( feature_alignment_requant , nrt = 3 , peptop = 3 )
#' dim(msexp)
#' sum(rownames(msexp$Intensity) == rownames(msexp$pepinfo))
#' y = aggregatePeptides(msexp)
#' sum(rownames(y$Intensity) == rownames(y$pepinfo))
#' dim(y)
#' stopifnot(table(table(y$pepinfo$transition_group_id)) == 247)
#' @seealso aggregateProtein
aggregatePeptides=function(msexp, FUN = sum)
{
  # fixing incopatibility with data.table
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  tomatrix = function(x){
    rownames(x) = x[,1]
    x = x[,-1]
    return(as.matrix(x))
  }
  #check sorting
  stopifnot(rownames(msexp$Intensity)==rownames(msexp$pepinfo))
  
  aggval = list( msexp$pepinfo$transition_group_id)
  intensity = aggregate(msexp$Intensity,by = aggval,FUN=FUN)
  res = msexp
  res$Intensity = tomatrix(intensity)
  x = aggregate(msexp$mz,by=aggval , function(x){x[1]})
  res$mz = tomatrix(x)
  x = aggregate(msexp$rt,by=aggval , function(x){x[1]})
  res$rt = tomatrix(x)
  x = aggregate(msexp$score,by=aggval , function(x){x[1]})
  res$score = tomatrix(x)
  
  res$pepinfo = msexp$pepinfo[,-match("aggr_Fragment_Annotation",colnames(msexp$pepinfo))]
  tmp = duplicated(res$pepinfo) 
  res$pepinfo =res$pepinfo[!tmp,]
  stopifnot(res$pepinfo$transition_group_id == rownames(res$Intensity))
  rownames(res$pepinfo) = res$pepinfo$transition_group_id
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(res)
}
#' aggregate peptides - given msexperiment with transitons
#' @export 
#' @param msexp - msexperiment
#' @param FUN aggregation function to use 
#' @examples
#' data( feature_alignment_requant )
#' msexp = loadTransitonsMSExperiment( feature_alignment_requant , nrt = 3 , peptop = 3 )
#' x = msexp
#' x = removeDecoys( x )
#' table( table( x$pepinfo$ProteinName ) )
#' y = aggregateProteins( x )
#' stopifnot(rownames(y$pepinfo$ProteinName) == rownames(y$Intensity))
#' stopifnot( table( table( y$pepinfo$ProteinName ) ) == 50 )
#' 
#' @seealso aggregatePeptide
aggregateProteins=function( msexp, FUN = sum )
{
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  
  tomatrix = function(x){
    rownames(x) = x[,1]
    x = x[,-1]
    return(as.matrix(x))
  }
  res = msexp
  stopifnot(rownames(res$Intensity)==rownames(res$pepinfo))
  
  aggval = list( res$pepinfo$ProteinName)
  intensity = aggregate(res$Intensity,by= aggval,FUN=FUN)
  res$Intensity = tomatrix(intensity)
  #x = aggregate(msexp$mz,by=aggval,function(x){x[1]})
  #set null since there are no peptide values anymore
  res$mz = NULL
  res$rt = NULL
  # compute mean score per protein (TODO: - think about requant)
  x = aggregate(res$score,by=aggval,mean)
  res$score = tomatrix(x)
  res$pepinfo = res$pepinfo[, c("ProteinName", "decoy")]
  res$pepinfo = res$pepinfo[order(res$pepinfo$ProteinName),]

  tmp = duplicated( res$pepinfo$ProteinName ) 
  res$pepinfo = res$pepinfo[ !tmp , ]
  
  stopifnot(res$pepinfo$ProteinName == rownames(res$Intensity))
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(res)
}
