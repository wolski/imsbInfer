#' plot log2(x/y) against retention time
#' @param data data.frame with intensities
#' @param RT rt information (same langth as nrow(data))
#' @param main main title
#' @param ylim y image extend
#' @export
#' @examples
#' data(SDat)
#' rtord = orderByRT(SDat)
#' pairsRatio( rtord$Intensity, rtord$RT )
#' @seealso \code{\link{pairs}} and \code{\link{pairsDifference}}
pairsRatio = function(data,RT,main="",ylim=c(-2,2)){
  pairs(data, panel = function(x,y,...){
    y = log2(x / y)
    r <- points(RT[!is.na(y)] , y[!is.na(y)],pch=".",... )
    abline(h = 0,col=2)
  }
  ,xlim=range(RT)
  ,ylim=ylim
  , lower.panel=NULL
  ,main = main
  )
}
#' plot x-y against retention time
#' @export
#' @param dataf dataframe
#' @param RT retention time vector
#' @param main title
#' @param ylim axis extend
#' @param maxPanel maximum number of datasets to display. dim(dataf)[2] > maxPanel, maxPanel columns are sampled
#' @examples
#' data(SDat)
#' rtord = orderByRT(SDat)
#' pairsDifference( asinh(rtord$Intensity), rtord$RT , ylim=NULL)
#' @seealso \code{\link{pairs}} and \code{\link{pairsRatio}}
pairsDifference <- function(dataf, RT, main="", ylim=NULL, maxPanel=8 ){
  stopifnot(dim(dataf)[2] > 1)
  if(dim(dataf)[2] > maxPanel){
    idx = sample(1:dim(dataf)[2],maxPanel)
    dataf = dataf[,idx]
  }
  if(is.null(ylim)){
    nr = dim(dataf)[2]
    tmp = NULL
    for(i in 1:nr){
      for(j in 1:nr){
        if(i != j){
          tmp = range(c(tmp,(dataf[,i] - dataf[,j])), na.rm=TRUE)
        }
      }
    }
    ylim = tmp
  }
  print(ylim)
  
  pairs(dataf, panel = function(x,y,...){
    y = x - y
    r <- points(RT[!is.na(y)] , y[!is.na(y)],pch=".",... )
    abline(h = 0,col=2)
  }
  ,xlim=range(RT)
  ,ylim=ylim
  ,lower.panel=NULL
  ,main = main
  )
}




