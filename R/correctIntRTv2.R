#' running robust scaling of arefw
#' @export
#' @examples
#' res = c(rnorm(1000,0,1),rnorm(2000,0,3))
#' res2 = runrobscale(res)
#' plot(res,type="p",pch="x",col=1)
#' lines(res2$runmed,col=3)
#' points(res2$scaled, pch=".",cex=3,col=2)
#' @seealso correctIntRTv2 for context
runrobscale = function(arefw,k=101,scale=TRUE){
  medianref = runmed( arefw , k=k ,endrule="constant")
  # adjust intensities
  # center around 0
  scalefactor = medianref
  mediana1w = arefw - scalefactor
  if(scale){
    madref <- runFun(arefw,k=k)
    mediana1w <- mediana1w / madref
  }
  return(list("scaled"  = mediana1w, "runmed" = medianref))
}
#' running total ion count scaling (TIC)
#' @export
#' @examples
#' res = c(rnorm(1000,3,2),rnorm(2000,1,1))
#' res2 = runTIC(res)
#' plot(res,type="p",pch="x",col=1)
#' lines(1:length(res),res2$tic,col=3)
#' points(res2$scaled, pch=".",cex=3,col=2)
#' length(res2$scaled)
#' 
#' @seealso correctIntRTv2 for context
runTIC = function(arefw,k=101,scale=TRUE){
  madref <- runFun(arefw,k=k,function(x,...){return(mean(x))})
  mediana1w <- arefw / madref
  return(list("scaled"  = mediana1w, "tic" = madref))
}

#' @export
correctIntRTv2 <- function(obj, ... ){
  UseMethod('correctIntRTv2')
}
#' correct Intensity RT using func
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
#' @examples
#' res = c(rnorm(1000,0,1),rnorm(2000,0,3))
#' res[sample(1:length(res),100)] = NA
#' rto = as.numeric(1:length(res))
#' res2 = correctIntRTv2(res ,rto )
correctIntRTv2.default = function(obj, rto , scale=FALSE, plot=TRUE, k=501, func=runrobscale, ...){
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
  bb[!idxref] = resScale$runmed
  return(bb)
}
#' correct intensity over RT for entire msexperiment
#' 
#' @note reorders all entries in experiment according to RT, finds sample with fewest NA's if there are many picks that one
#' with largest median as reference.
#' @param experiment - object of class msexperiment
#' @param k - smoothing with
#' @return msexperiment object with RT normalized intensities
#' @export
#' @author Witold Wolski \email{wolski@@gmail.com}
#' @examples
#' data(SDat)
#' res = correctIntRTv2(SDat)
correctIntRTv2.msexperiment = function(obj , k=501,plot=FALSE , scale=FALSE, func = runrobscale, ... )
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

