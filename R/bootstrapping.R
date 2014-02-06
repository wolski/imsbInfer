#' bootstrapping1 
#' @param nr number of samplings
#' @return 3 D matrix where dim()[1] == nr
#' @export
computeBoost <- function(pln,nr=250)
{
  #result1=ffunc(pln)
  bla <- dim(pln)
  dd <- array(0,dim=c(nr,bla[1],bla[2]))
  for(i in 1:nr){
    print(i)
    ps <- pln[sample(dim(pln)[1],replace=TRUE),]
    dd[i,,] <- ffunc(as.matrix(ps))
  }
  return(dd)
}
#' bootstrapping 2 
#' @param pln data matrix
#' @param g1 group 1 indices
#' @param g2 group 2 indices
#' @export
computeBoost2 <- function(pln,g1,g2,ffunc,nr=250)
{
  result1 <- ffunc(pln)
  bla <- dim(result1)
  dd <- array(0,dim=c(nr,bla[1],bla[2]))
  for(i in 1:nr){
    print(i)
    ps <- pln[c(
      g1[sample(length(g1),replace=TRUE)]
      ,g2[sample(length(g2),replace=TRUE)]
    ),]
    dd[i,,] <- ffunc(as.matrix(ps))
  }
  return(dd)
}
#' Methods for selecting significant 
#'
#' works in conjunction with the 2 bootstrap methods
selectSignificant2 <- function(dd, siglevel = 0.7, minPcor=0)
{
  tmpq5 <- apply(dd,c(2,3),quantile,1-siglevel)
  tmpq95 <- apply(dd,c(2,3),quantile,siglevel)
  diag(tmpq5) <- 0
  diag(tmpq95) <- 0
  
  indic <- (tmpq5 > minPcor | tmpq95 < -minPcor)
  
  tmpmedian <- apply(dd,c(2,3),median)
  tmpmedian[!indic] <- 0
  diag(tmpmedian)=0
  return(tmpmedian)
}
#' Methods for selecting significant partial correlation
#'
#' Does statistics over bootstrap statistics
#' @export
selectSignificant <- function(dd, siglevel = 0.8, minPcor=0)
{
  tmpq5 <- apply(dd,c(2,3),quantile,1-siglevel)
  tmpq95 <- apply(dd,c(2,3),quantile,siglevel)
  diag(tmpq5) <- 0
  diag(tmpq95) <- 0
  
  indic <- (tmpq5 > 0 | tmpq95 < 0)
  
  tmpmedian <- apply(dd,c(2,3),median)
  tmpmedian[!indic] <- 0
  diag(tmpmedian)=0
  indic <- (tmpmedian > minPcor | tmpmedian < -minPcor)
  tmpmedian[!indic] <- 0
  return(tmpmedian)
}
