#' correct intensity over RT for entire msexperiment
#'
#' @export
correctIntRTv1 <- function(obj, ... ){
  UseMethod('correctIntRTv1')
}
#' correct intensity over RT for entire msexperiment
#' @param obj  reference data - ordered by retention time
#' @param data to correct
#' @param rto  retention time (ordered)
#' @param k  smoothing with
#' @param FUN running function - for equal median use runmed, for equal max use i.e. runFUN(x,k=300,max)
#' @return corrected data
#' @export
#' @author Witold Wolski \email{wolski@@gmail.com}
#' @examples
#' data(SDat)
#' experiment = removeDecoys(SDat)
#' experiment = orderByRT(experiment)
#' cor = correctIntRTv1(experiment$Intensity[,1], experiment$Intensity[,2],experiment$RT,k=51)
#' cor = correctIntRTv1(experiment$Intensity[,1], experiment$Intensity[,2],experiment$RT,k=51,FUN=function(x,k=k){print(k);runFun(x,k=k,mean)})
correctIntRTv1.default <- function(obj, data, rto , plot=TRUE,k=501, FUN = function(x,k=k){runmed(x,k=k,endrule="constant")}, ...){
  aref = obj  
  a1=data
  #remove missing values
  idxref = is.na(aref) |is.infinite(aref)
  idxs = is.na(a1) | is.infinite(a1) 
  idx =  !(idxref | idxs)
  a1w = a1[idx]
  arefw = aref[idx]
  rtow = rto[idx]
  
  #diffa1ref  = a1w - arefw
  medianref = FUN( arefw , k=k)# ,endrule="constant")
  mediana1w = FUN( a1w ,k=k)# ,endrule="constant")
  # adjust intensities
  scalefactor = medianref/mediana1w
  a1wc = a1w * scalefactor
  
  if(plot){
    plot(rtow , arefw , pch=".", cex=0.4 , col="gray" )
    points( rtow , a1w ,pch=".",cex=0.4,col="gray")
    lines(rtow , medianref,col="red",lwd=4)
    lines(rtow , mediana1w,col="blue",lwd=2)
    #abline(h=0,col=2)
    mediana1wc = FUN( a1wc ,k = k )#endrule="constant")
    lines(rtow,mediana1wc,col="black",lwd=2)
    legend("topleft",legend=c("ref","before","after"),col=c("red","blue","black"),lty=c(1,1,1))
  }
  bb  = rep(NA, length(idxs))
  bb[idx] = a1wc
  return(bb)
}
#' correct intensity over RT for entire msexperiment
#'
#' @export
# select reference dataset
correctIntRTv1.matrix = function(obj,rt,k=501,plot=F,...)
{
  intensity = obj
  nas = apply( intensity, 2 , function(x){sum(is.na(x))} )
  idx = which(nas == min(nas))
  # if more than one with few NA's than choose dataset with max median
  if(length(idx) > 1){
    ma = apply(intensity[,idx],2,median,na.rm=TRUE)
    id <-which(ma == max(ma))
    idx <- idx[id]
  }
  
  reference=intensity[,idx]
  for(i in 1:dim(intensity)[2])
  {
    intensV = intensity[,i]
    corrected = correctIntRTv1( obj  = unlist(reference), data = unlist(intensV) , rto = rt , plot=plot , k=k )
    intensity[,i] = corrected
  }
  return(intensity)
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
#' SDatal = correctIntRTv1(SDat,plot=F)
#' 
correctIntRTv1.msexperiment = function(obj,k=501,plot=F, ...){
  experiment = obj
  experiment = removeDecoys(experiment)
  experiment = orderByRT(experiment)
  experiment$Intensity = correctIntRTv1(obj = experiment$Intensity,rt = experiment$RT,k=501,plot=F,...)
  return(experiment)
}

