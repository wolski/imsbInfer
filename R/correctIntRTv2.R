#' @export
correctIntRTv2 <- function(obj, ... ){
  UseMethod('correctIntRTv2')
}
#' correct Intensity RT
#' Computes difference between aref and data and adjust data
#'  that the median of differences is equal 0
#' @aliases correctIntRTv2
#' @param aref  reference data - ordered by retention time
#' @param data to correct
#' @param rto  retention time
#' @param k  smoothing with
#' @return corrected data
#' @export
#' @author Witold Wolski \email{wolski@@gmail.com}
correctIntRTv2.default = function(aref, rto , scale=FALSE, plot=TRUE, k=501){
  idxref = is.na(aref) | is.infinite(aref)
  arefw = aref[!idxref]
  rtow = rto[!idxref]
  length(arefw)
  #diffa1ref  = a1w - arefw
  medianref = runmed( arefw , k=k ,endrule="constant")
  # adjust intensities
  scalefactor = medianref
  mediana1w = arefw - scalefactor
  if(scale){
    madref <- runMAD(arefw,k=k)
    mediana1w <- mediana1w / madref
  }
  
  if(plot){
    par(mfrow=c(2,1))
    plot(rtow , arefw , pch=".", cex=0.4 , col="gray" ,ylim=c(-10,15))
    lines(rtow , medianref,col="red",lwd=2)
    plot(rtow , mediana1w,col="blue",pch=".")
    medianref2 = runmed( mediana1w , k=k ,endrule="constant")
    lines(rtow , medianref2,col="red",lwd=2)
  }
  bb  = rep(NA, length(idxref))
  bb[!idxref] = mediana1w
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
#' 
correctIntRTv2.msexperiment = function(experiment , k=501,plot=FALSE , scale=FALSE)
{
  experiment = removeDecoys(experiment)
  experiment = orderByRT(experiment)
  
  for(i in 1:dim(experiment$intensity)[2])
  {
    intensV = experiment$intensity[,i]
    corrected = correctIntRTv2( unlist(intensV) , (experiment$RT) , plot=plot , k=k, scale=scale )
    experiment$intensity[,i] = corrected
  }
  return(experiment)
}

