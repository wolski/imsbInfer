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
#' 
#' @export
split2table <- function(names,split="\\||\\_")
{
  cnamessplit <- strsplit(names,split)
  protnam <- NULL
  for(i in cnamessplit)
  {
    protnam <- rbind(protnam,i)
  }
  return(protnam)
}
#' getupper triangle from matrix 
#' @export
uppertriang <- function(mat){
  res<-mat[upper.tri(mat,diag=FALSE)] 
  return( c(unlist(res)) )
}
