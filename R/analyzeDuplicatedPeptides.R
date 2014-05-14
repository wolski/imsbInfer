#' analyse duplicated peptides
#'
#' @export
analyzeDuplicatedPeptides = function(dups,data){
  res = NULL
  count = 1
  nameshead = c("p1","p2","cor","rt1","rt2","medianRTDiff","madDiffRT")
  res = matrix(NA ,nrow = length(dups)*3 ,ncol=length(nameshead) )
  dim(res)
  
  for(dup in dups){
    #dup = dups[1]
    duplicated <- which(data$pepinfo$sequence %in% dup)
    ld =length(duplicated)
    print(dup)
    print((ld))
    for(i in 1:ld){
      for(j in i:ld){
        if(i!=j){
          tmp=cor(t(data$intensity[duplicated[c(i,j)],]))
          cors = tmp[upper.tri(tmp)]
          x<-which(upper.tri(tmp),arr.ind=T)
          nams <- rownames(tmp)[x]
          medianrt = apply(data$rt,1,median)
          rowstoget = which(rownames(data$rt) %in% nams)
          RTDiff = data$rt[rowstoget[1],] - data$rt[rowstoget[2],]
          tmp=c(t(nams),cors,t( medianrt[rowstoget]), median(RTDiff), mad(RTDiff) )
          #names(tmp) = nameshead
          res[count,] = tmp
          count = count + 1
        }
      }
    }
  }
  res2= as.data.frame(res[1:count,],stringsAsFactors=FALSE)
  colnames(res2) = nameshead
  res2$cor=as.numeric(res2$cor)
  res2$rt1=as.numeric(res2$rt1)
  res2$rt2=as.numeric(res2$rt2)
  res2$madDiffRT=as.numeric(res2$madDiffRT)
  return(res2)
}
