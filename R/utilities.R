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
rmNarows <- function(plev, thresh=0 )
{
  x <- apply(plev,1,function(x){sum(is.na(x))})
  plev <- plev[-which(x>thresh),]
}
#' splits names and binds into matrix
#' @examples
#' dat = c("bla_ra0/2_run0","bla_ra1/2_run0","bla_ra2/2_run0")
#' split2table(dat,split="\\_|\\/")
#' @export
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
#' getupper triangle from matrix 
#' @export
uppertriang <- function(mat){
  res<-mat[upper.tri(mat,diag=FALSE)] 
  return( c(unlist(res)) )
}
#' running median absolute deviation (not quite efficient)
#'
#' @examples
#' x = rnorm(5000)
#' y = runMAD(x,k=501)
#' @export
#' @seealso  \code{\link{runmed}}
#' 
runMAD <- function(aref,k=301){
  m = k %/% 2
  N = length(aref)
  res = rep(0,N)
  for(i in 1:(N-k+1)){
    sub = aref[i:(i+k)]
    res[i + m]=mad(sub,na.rm=T)
  }  
  res[1:m] = rep(res[m+1],m)
  res[(N-m + 1):N] = rep(res[N-m],m)
  return(res)
}
