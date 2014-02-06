#' correct Intensity RT
#' Computes difference between aref and data and adjust data
#'  that the median of differences is equal 0
#' @param aref - reference data - ordered by retention time
#' @param data to correct
#' @param rto - retention time
#' @return corrected data
#' @export
#' @author Witold Wolski \email{wolski@@gmail.com}
correctIntRT = function(aref, data, rto , plot=TRUE){
  a1=data
  #remove missing values
  idxref = is.na(aref) |is.infinite(aref)
  idxs = is.na(a1) | is.infinite(a1) 
  idx =  !(idxref | idxs)
  a1w = a1[idx]
  arefw = aref[idx]
  rtow = rto[idx]
  
  diffa1ref  = a1w - arefw
  mediandiff = runmed( diffa1ref , k=501 )
  ac1w = ( a1w - mediandiff )
  if(plot){
    plot(rtow , diffa1ref , pch="." , cex=0.4 , col=1 )
    points( rtow , ac1w - arefw ,pch=".",cex=0.4,col=4)
    lines(rtow,mediandiff,col=3,lwd=2)
    abline(h=0,col=2)
  }
  bb  = rep(NA, length(idxs))
  #introduce back missing values
  bb[!idxs] = ac1w
  return(bb)
}
