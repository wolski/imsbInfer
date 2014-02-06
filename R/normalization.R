#' Normalization by total sum
#' The bla an the blo
NormalizeTotalSumOfPQMatrix = function(ProtQuantMatrix){
  maxSum <- max(colSums(ProtQuantMatrix))
  ScaleFactors <- 1/(apply(na.omit(ProtQuantMatrix),2,sum)/maxSum)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactors,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactors))
}
#' Normalization by total sum
#' The bla an the blo
#' @export
NormalizeWithMedianPQMatrix = function(ProtQuantMatrix){
  maxMedian <- max(apply(na.omit(ProtQuantMatrix),2,median))
  ScaleFactors <- 1/(apply(na.omit(ProtQuantMatrix),2,median)/maxMedian)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactors,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactors ))
  
}
#' Normalization by total sum
#' The bla an the blo
#' @export 
NormalizePQMatrixWithHouseKeepingProtein = function(ProtQuantMatrix, hk){
  MaxValueOfHK <- max(ProtQuantMatrix[hk,])
  ScaleFactorsHK <- 1/(ProtQuantMatrix[hk,]/MaxValueOfHK)
  nPQmatrix <- sweep(ProtQuantMatrix,2,ScaleFactors,"*")
  return(list(matrix=nPQmatrix, scalefactors=ScaleFactors ))
}
