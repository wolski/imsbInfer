#' pvalues for correlation matrix - Benjamin Hochberg - adjusted
#' @param x data matrix 
#' @param alternative (see cor.test)
#' @param method (see cor.test) 
#' @export
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
  return(list(res = res2, cor = corr))
}




