#' Normalization by total sum
#' @export
NormalizeTotalSumOfPQMatrix = function(ProtQuantMatrix){
  maxSum <- max(colSums(ProtQuantMatrix))
  ScaleFactors <- 1/(apply(na.omit(ProtQuantMatrix),2,sum)/maxSum)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactors,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactors))
}
#' Normalization by median
#' @export
NormalizeWithMedianPQMatrix = function(datamatrix){
  meds=apply( datamatrix , 2 , median , na.rm=TRUE )
  params = max( meds )/meds
  mednorm = sweep(tmp,2,params,"*")
  return(list(matrix=mednorm, scalefactors=params ))
  
}
#' Normalization by total sum
#' @export 
NormalizePQMatrixWithHouseKeepingProtein = function(ProtQuantMatrix, hk){
  MaxValueOfHK <- max(ProtQuantMatrix[hk,])
  ScaleFactorsHK <- 1/(ProtQuantMatrix[hk,]/MaxValueOfHK)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactorsHK,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactorsHK ))
}
