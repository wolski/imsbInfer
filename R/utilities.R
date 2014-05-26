#' Depracated 
#' @export
getMolfun <- function(pfilename){
  pgl <- read.csv(pfilename , stringsAsFactors=FALSE)
  gofunction <- pgl$GOMolecularFunction
  dd <- lapply(gofunction , function(x){d <- strsplit(x,";"); return(d[[1]][1])})
  dd <- lapply(dd , function(x){d <- strsplit(x,","); return(d[[1]][1])})
  tmp <- (unlist(dd))
  return(list(ids=pgl$Mapped.IDs, molfun=tmp))
}
#' Removes rows with nas from matrix
#' @export
rmNArows <- function(plev, thresh=0 )
{
  x <- apply(plev,1,function(x){sum(is.na(x))})
  plev <- plev[-which(x>thresh),]
}
#' splits names and binds into matrix
#' @export
#' @examples
#' dat = c("bla_ra0/2_run0","bla_ra1/2_run0","bla_ra2/2_run0")
#' split2table(dat,split="\\_|\\/")
split2table <- function(names,split="\\||\\_")
{
  cnamessplit <- strsplit(names,split)
  protnam <- matrix(NA, ncol=length(cnamessplit[[1]]),nrow=length(cnamessplit))
  print(dim(protnam))
  for(i in 1:length(cnamessplit))
  {
    protnam[i,] <- cnamessplit[[i]]
  }
  return(protnam)
}
#' get values of upper triangle from matrix 
#' @export
#' @examples
#' t = matrix(1:25,ncol=5)
#' uppertriang(t)
uppertriang <- function(mat){
  res<-mat[upper.tri(mat,diag=FALSE)] 
  return( c(unlist(res)) )
}
#' running function (default median absolute deviation) 
#' @param func default med but can be any function taking a vector and returning a summary
#' @examples
#' x = rnorm(500)
#' y = runFun(x,k=51)
#' hist(y)#[500:490]
#' y2 = runFun(x,k=51,func=median)
#' plot(x,pch="*")
#' lines(y2,col=2,lwd=3)
#' lines(-y,col=3,lwd=3)
#' tic = runFun(x,k=51,func=function(x,...){mean(x)})
#' plot(x,pch=".")
#' abline(h=0,col=2)
#' lines(tic,col=3,lwd=3)
#' @export
#' @seealso  \code{\link{runmed}}
runFun <- function(aref,k=301,func=mad)
  {
  m = k %/% 2
  N = length(aref)
  res = rep(NA,N)
  for(i in 1:(N-k)){
    sub = aref[i:(i+k)]
    res[i + m]=func(sub,na.rm=TRUE)
  }
  res[1:m] = rep(res[m+1],m)
  res[(N-m ):N] = rep(res[N-m-1],m+1)
  return(res)
}
