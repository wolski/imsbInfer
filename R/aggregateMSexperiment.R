
#' aggregate peptides - given msexperiment with transitons
#' @export
#' @examples
#' NULL
#' @seealso aggregateProtein
#' 
aggregatePeptides=function(msexp, FUN = sum)
{
  tomatrix = function(x){
    rownames(x) = x[,1]
    x = x[,-1]
    return(as.matrix(x))
  }
  
  #check sorting
  #sum(rownames(msexp$Intensity)==rownames(msexp$pepinfo))
  
  aggval = list( msexp$pepinfo$transition_group_id)
  intensity = aggregate(msexp$Intensity,by= aggval,FUN=FUN)

  msexp$Intensity = tomatrix(intensity)
  x = aggregate(msexp$mz,by=aggval,function(x){x[1]})
  msexp$mz = tomatrix(x)
  x = aggregate(msexp$rt,by=aggval,function(x){x[1]})
  msexp$rt = tomatrix(x)
  x = aggregate(msexp$score,by=aggval,function(x){x[1]})
  msexp$score = tomatrix(x)
  
  msexp$pepinfo = msexp$pepinfo[,-match("aggr_Fragment_Annotation",colnames(msexp$pepinfo))]
  tmp = duplicated(msexp$pepinfo) 
  msexp$pepinfo =msexp$pepinfo[!tmp,]
  return(msexp)
}


#' aggregate peptides - given msexperiment with transitons
#' @export 
#' @examples
#' NULL
#' @seealso aggregatePeptide
aggregateProteins=function(msexp, FUN = sum)
{
  tomatrix = function(x){
    rownames(x) = x[,1]
    x = x[,-1]
    return(as.matrix(x))
  }
  sum(rownames(msexp$Intensity)==rownames(msexp$pepinfo))
  aggval = list( msexp$pepinfo$ProteinName)
  intensity = aggregate(msexp$Intensity,by= aggval,FUN=FUN)
  msexp$Intensity = tomatrix(intensity)
  #x = aggregate(msexp$mz,by=aggval,function(x){x[1]})
  msexp$mz = NULL
  msexp$rt = NULL
  x = aggregate(msexp$score,by=aggval,mean)
  msexp$score = tomatrix(x)
  
  msexp$pepinfo = msexp$pepinfo[,-match(c("aggr_Fragment_Annotation","transition_group_id") , colnames( msexp$pepinfo ) )]
  tmp = duplicated(msexp$pepinfo) 
  msexp$pepinfo =msexp$pepinfo[!tmp,]
  return(msexp)
}
