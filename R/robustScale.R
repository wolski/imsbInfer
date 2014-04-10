#' robust scaling
#' uses median an mad instead of mean and row
#' applies the scaling to the columns (samples) by default
#' @export
robustscale <- function(data, dim=2, center=TRUE, scale=TRUE){
  if(center){
    medians <- apply(data,dim,median,na.rm=TRUE)
    data = sweep(data,dim,medians,"-")
  }
  if(scale){
    mads <- apply(data,dim, mad,na.rm =TRUE)
    data = (sweep(data,dim,mads,"/"))
  }
  return(data)
}
