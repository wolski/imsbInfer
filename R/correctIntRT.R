#' correctIntRT
#'
#' @export
correctIntRT <- function(obj, ...){
  UseMethod('correctIntRT')
}
#' correct Intensity RT
#' Computes difference between aref and data and adjust data
#'  that the median of differences is equal 0
#' @param aref  reference data - ordered by retention time
#' @param data to correct
#' @param rto  retention time
#' @param k  smoothing with
#' @return corrected data
#' @export
#' @S3method correctIntRT default
#' @author Witold Wolski \email{wolski@@gmail.com}
correctIntRT.default = function(aref, data, rto , plot=TRUE,k=501){
  a1=data
  #remove missing values
  idxref = is.na(aref) |is.infinite(aref)
  idxs = is.na(a1) | is.infinite(a1) 
  idx =  !(idxref | idxs)
  a1w = a1[idx]
  arefw = aref[idx]
  rtow = rto[idx]
  
  #diffa1ref  = a1w - arefw
  medianref = runmed( arefw , k=k ,endrule="constant")
  mediana1w = runmed( a1w ,k = k,endrule="constant")
  # adjust intensities
  scalefactor = medianref/mediana1w
  a1wc = a1w * scalefactor
  
  if(plot){
    plot(rtow , arefw , pch=".", cex=0.4 , col="gray" )
    points( rtow , a1w ,pch=".",cex=0.4,col="gray")
    lines(rtow , medianref,col="red",lwd=2)
    lines(rtow , mediana1w,col="blue",lwd=2)
    abline(h=0,col=2)
    mediana1wc = runmed( a1wc ,k = k,endrule="constant")
    lines(rtow,mediana1wc,col="black",lwd=2)
    legend("topleft",legend=c("ref","before","after"),col=c("red","blue","black"),lty=c(1,1,1))
  }
  bb  = rep(NA, length(idxs))
  bb[!idxs] = a1wc
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
#' @S3method correctIntRT msexperiment
#' @author Witold Wolski \email{wolski@@gmail.com}
#' 
correctIntRT.msexperiment = function(experiment,k=501){
  experiment = removeDecoys(experiment)
  dim(experiment)
  
  #as reference choose run with fewest NA's
  firstDataColumn = 2

  # select reference dataset
  nas = apply( experiment$intensity, 2 , function(x){sum(is.na(x))} )
  idx = which(nas == min(nas))
  # if more than one with few NA's than choose dataset with max median
  if(length(idx) > 1){
    ma = apply(experiment$intensity[,idx],2,median,na.rm=TRUE)
    id <-which(ma == max(ma))
    idx <- idx[id]
  }
  reference=intens[,idx]
  
  #compute median retention time
  RT = apply(experiment$rt,1,median)
  #order relevant data by retention time
  rto = order(RT)
  reference = reference[rto]
  RT = RT[rto]
  experiment$intensity = experiment$intensity[rto,]
  experiment$score = experiment$score[rto,]
  experiment$rt = experiment$rt[rto,]
  
  for(i in 1:experiment$intensity)
  {
    intensV = intens[[i]]
    corrected = correctIntRT( unlist(reference), unlist(intensV) , (RT) , plot=F , k=k )
    intens[[i]] = corrected
  }
  experiment$intensity = intens
}
