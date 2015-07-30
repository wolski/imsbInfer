#' filter rows by qscore
#' @description require that at least n features have a mscore less than maxscore
#' @param msexp object of class msexperiment
#' @param scorethresh - score threshold
#' @param n - number of samples which need to have a score for this feature below the threshold
#' @return msexperiment 
#' @export
#' @examples
#' NULL
filterByScore = function(msexp, scorethresh = 0.001, n = 2)
{
  tmp = apply(msexp$score,1,function(x){sort(x)[n]})
  tmp[is.na(tmp)] = 2
  return(subset(msexp, tmp < scorethresh))
}
