
#' @export
correctIntRTv2 <- function(obj, ... ){
  UseMethod('correctIntRTv2')
}
#' correct Intensity RT using func, i.e. runrobscale (default)
#' @aliases correctIntRTv2
#' @param obj to correct
#' @param rto  retention time
#' @param plot show diagnostic plot
#' @param scale should scaling be applied
#' @param k smoothing with (see runmed)
#' @param func function to apply for scaling see runrobscale
#' @return corrected data
#' @export
#' @author Witold Wolski 
#' @seealso runrobscale runTICscale
#' @examples
#' res = c(rnorm(1000,0,1),rnorm(2000,0,3))
#' res[sample(1:length(res),100)] = NA
#' rto = as.numeric(1:length(res))
#' res2 = correctIntRTv2(res ,rto )
correctIntRTv2.default = function(obj, rto ,k=501, scale=FALSE, plot=TRUE,  func=runrobscale, ...){
  aref= obj
  idxref = is.na(aref) | is.infinite(aref)
  arefw = aref[!idxref]
  rtow = rto[!idxref]
  
  resScale = func(arefw,k=k,scale=scale)
  
  if(plot){
    par(mfrow=c(2,1))
    plot(rtow , arefw , pch=".", cex=0.4 , col="gray" ,ylim=c(-10,15))
    lines(rtow , resScale$runmed,col="red",lwd=2)
    points(rtow , resScale$scaled,col="blue",pch=".")
  }
  bb  = rep(NA, length(idxref))
  bb[!idxref] = resScale$scaled
  return(bb)
}
#' correct intensity over RT for entire msexperiment
#' 
#' @note reorders all entries in experiment according to RT, finds sample with fewest NA's if there are many picks that one
#' with largest median as reference.
#' @param experiment - object of class msexperiment
#' @param k - smoothing with
#' @param scale  should scaling also be applied
#' @param fun function to use fore scaling: \code{\link{runrobscale}} or \code{\link{runTICscale}}
#' @return msexperiment object with RT normalized intensities
#' @export
#' @author Witold Wolski \email{wolski@@gmail.com}
#' @examples
#' data(SDat)
#' res = correctIntRTv2(SDat)
correctIntRTv2.msexperiment = function(obj , k=501 ,  scale=FALSE, func = runrobscale,plot=FALSE , ... )
{
  experiment = obj
  experiment = removeDecoys(experiment)
  experiment = orderByRT(experiment)
  
  for(i in 1:dim(experiment$Intensity)[2])
  {
    intensV = experiment$Intensity[,i]
    corrected = correctIntRTv2( unlist(intensV) , (experiment$RT) , plot=plot , k=k, scale=scale, func = runrobscale )
    experiment$Intensity[,i] = corrected
  }
  return(experiment)
}

