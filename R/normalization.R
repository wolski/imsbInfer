#' Normalization by total sum
#' @param ProtQuantMatrix matrix of intensities
#' @export
normalizeTotalSum = function(ProtQuantMatrix){
  maxSum <- max(colSums(ProtQuantMatrix))
  ScaleFactors <- 1/(apply(na.omit(ProtQuantMatrix),2,sum)/maxSum)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactors,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactors))
}
#' Normalization by median
#' @param datamatrix matrix of intensities
#' @export
normalizeWithMedian = function(datamatrix){
  meds=apply( datamatrix , 2 , median , na.rm=TRUE )
  params = max( meds )/meds
  mednorm = sweep(tmp,2,params,"*")
  return(list(matrix=mednorm, scalefactors=params ))
}
#' Normalization by total sum
#' @param ProtQuantMatrix matrix of intensities
#' @param hk indices of house keeping proteins
#' @export 
normalizeWithHouseKeepingProtein = function(ProtQuantMatrix, hk){
  MaxValueOfHK <- max(ProtQuantMatrix[hk,])
  ScaleFactorsHK <- 1/(ProtQuantMatrix[hk,]/MaxValueOfHK)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactorsHK,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactorsHK ))
}
