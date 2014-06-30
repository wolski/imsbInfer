#' Removes rows with NA's from matrix
#' @export
#' @return matrix
#' @examples
#' x = matrix(rnorm(10*10),ncol=10)
#' dim(x)
#' x[3,3] = NA
#' x = rmNArows(x)
#' dim(x)
rmNArows <- function(plev, thresh=0 )
{
  x <- apply(plev,1,function(x){sum(is.na(x))})
  plev <- plev[-which(x>thresh),]
}
#' splits names and creates a matrix
#' @export
#' @return matrix
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
#' x = c(x,rnorm(1000,3,2))
#' x = c(x,runif(1000,4,6))
#' y = runFun(x,k=51,func=mad)
#' hist(y)#[500:490]
#' y2 = runFun(x,k=51,func=median)
#' plot(x,pch="*")
#' lines(y2,col=2,lwd=3)
#' lines(y2+y,col=3,lwd=3)
#' lines(y2-y,col=3,lwd=3)
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
