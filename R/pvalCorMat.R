#' pvalues for correlation matrix - Benjamin Hochberg - adjusted
#' tests for significance of correlations among rows
#' transpose your matrix if you want to have it among columns
#' @param x data matrix 
#' @param alternative (see cor.test)
#' @param method (see cor.test) 
#' @export
#' @examples
#'  
#' mat = matrix(rnorm(10*20),ncol=10)
#' res = pvalCorMat(mat)
#' image(res$pval)
#' image(res$cor)
pvalCorMat = function(x, alternative = "two.sided", method = "spearman"){
  resP = matrix(0,nrow=dim(x)[1],ncol=dim(x)[1])
  corr = resP
  res2 = resP
  for(i in 1:dim(x)[1])
  {
    for(j in 1:dim(x)[1])
    {
      if(i < j)
      {
        test = cor.test(x[i,],x[j,],alternative,method)
        resP[i,j] = test$p.value
        corr[i,j] = test$estimate
      }
    }
  }
  ## adjust p values 
  res2[upper.tri(res2)] = (p.adjust(resP[upper.tri(resP)],method="BH"))
  res2[lower.tri(res2 , diag=TRUE)] = 1
  
  res2[lower.tri(res2)] = t(res2)[lower.tri(res2)]
  corr[lower.tri(corr)] = t(corr)[lower.tri(corr)]
  return(list(pval = res2, cor = corr))
}

#' dist with freely choosable distance function
#' @param x data
#' @param func function taking 2 arrays x, y
#' @param init how to initialize the output matrix
#' @param diag should the diagonal be also computed
#' @export
#' @examples
#' mat = matrix(rnorm(10*5000),ncol=10)
#' redist = distmy(mat,function(x,y){mean(abs(x-y))},init=NA,diag=FALSE)
#' image(redist)
#' redist = distmy(mat,cor,init=0,diag=FALSE)
#' image(redist)
#' redist = distmy(mat,function(x,y){ks.test(x,y)$p.value},init=1,diag=TRUE)
#' image(redist)
#' hist(uppertriang(redist))
#' range(uppertriang(redist))
#' which(redist < 0.05 , arr.ind = TRUE)
distmy = function( x, func, init=NA , diag = TRUE)
{
  f = NULL
  if(diag){
    f = function(i,j){i<=j}
  }else{
    f = function(i,j){i<j}
  }
  nout = dim(x)[2]
  resP = matrix(init,nrow=nout,ncol=nout)
  for(i in 1:nout){
    for(j in 1:nout){
      if(f(i,j)){
        test = func(x[,i],x[,j])
        resP[i,j] = test
      }
    }
  }
  resP[lower.tri(resP)] = t(resP)[lower.tri(resP)]
  rownames(resP) = colnames(x)
  colnames(resP) = colnames(x)
  return(resP)
}




