#' running robust scaling of arefw
#' @export
#' @examples
#' res = c(rnorm(1000,0,1),rnorm(2000,4,3))
#' res2 = runrobscale(res)
#' par(mfrow=c(2,1))
#' plot(res,type="p",pch="x",col=1,cex=0.3)
#' lines(res2$runmed,col=3)
#' 
#' y = runFun( res2$scaled, k=51, func=mad )
#' #hist(y)
#' y2 = runFun(res2$scaled,k=51,func=median)
#' plot(res2$scaled,pch="*")
#' lines(y2,col=2,lwd=3)
#' lines(y2+y,col=3,lwd=3)
#' lines(y2-y,col=3,lwd=3)
#' 
#' @seealso correctIntRTv2 for context
runrobscale = function(arefw,k=101,scale=TRUE){
  medianref = runmed( arefw , k=k ,endrule="constant")
  # adjust intensities
  # center around 0
  scalefactor = medianref
  mediana1w = arefw - scalefactor
  madref <- NULL
  if(scale){
    madref <- runFun(arefw,k=k)
    mediana1w <- mediana1w / madref
  }
  return(list("scaled"  = mediana1w, "runmed" = medianref, "runmad" = madref))
}
#' running total ion count scaling (TIC)
#' @export
#' @examples
#' res = c(rnorm(1000,3,2),rnorm(2000,8,1))
#' res2 = runTICscale(res)
#' plot(res,type="p",pch=".",col=1,cex=0.5)
#' lines(1:length(res),res2$tic,col=3)
#' points(res2$scaled, pch=".",cex=3,col=2)
#' 
#' @seealso correctIntRTv2 for context
runTICscale = function(arefw,k=101,scale=NULL){
  madref <- runFun(arefw,k=k,mean)
  mediana1w <- arefw / madref
  return(list("scaled"  = mediana1w, "tic" = madref))
}

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
#' @param fun function to use fore scaling: \code{\link{runrobscale}} or \code{\link{runTIC}}
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

