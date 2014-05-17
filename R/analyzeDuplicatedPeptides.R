#' get a list of random matches
#' 
#' @export
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$sequence
#' res = getUniquePeptides(rownames)
#' res[1:10]
#' sum(duplicated(res))
#' res2 = analyzeDuplicated( data , list(res[1:10] ))
#' dim(res2)
#' @seealso \code{\link{analyzeDuplicated}} for contex
getUniquePeptides = function(  rownames )
{
  dups = rownames[which(duplicated(rownames))]
  res = setdiff(rownames,dups)
  res = which( rownames %in% res)
  return(res)
}
#' randomly match peptides
#' @export
#' @examples
#' data(SDat)
#' res = randomPeptidePairs(1000,dim(SDat)[1])
#' res2 = analyzeDuplicated( data , res)
#' dim(res2)
#' @param number of pairs
#' @seealso \code{\link{analyzeDuplicated}}
randomPeptidePairs = function(n, dsize){
  res = list(n)
  tmp = 1:dsize
  for(i in 1:n){
    res[[i]] <- sample(tmp,2)
  }
}
#' get a list of duplicates
#' 
#' @export
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$sequence
#' res = getListOfMatches(rownames)
#' rownames[res[[1]]]
#' @seealso \code{\link{analyzeDuplicated}} for contex
getListOfMatches = function(  rownames ){
  dups = rownames[which(duplicated(rownames))]
  res = list(length(dups))
  for(i in 1:length(dups)){
    res[[i]]= which(rownames %in% dups[i])
  }
  return(res)
}
#' analyse duplicated peptides - same peptide different charge
#'
#' @export
#' @examples
#' data(SDat)
#' res = analyzeDuplicatedPeptides(SDat)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' #same for proteins
#' rownames = as.character(SDat$pepinfo$ProteinName)
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedProteins}}
analyzeDuplicatedPeptides = function(data){
  rownames = data$pepinfo$sequence
  dups=getListOfMatches( rownames )
  analyzeDuplicated( data ,dups )
}
#' analyse duplicated protein - same protein id different peptide
#'
#' @export
#' @examples
#' data(SDat)
#' res = analyzeDuplicatedProtein(SDat)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' #same for proteins
#' rownames = as.character(SDat$pepinfo$ProteinName)
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedPeptides}}
#' @seealso analyzeDuplicated
#' @param data
analyzeDuplicatedProteins = function(data){
  rownames = data$pepinfo$ProteinName
  dups=getListOfMatches( rownames )
  analyzeDuplicated( data ,dups )
}
#' analyse duplicated peptides/proteins/lines
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$sequence
#' dups = getListOfMatches(rownames)
#' res = analyzeDuplicated(SDat , dups[1:25])
#' dim(res)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' 
#' rownames = SDat$pepinfo$sequence
#' nondups = getUniquePeptides(rownames)
#' length(nondups)
#' res = analyzeDuplicated(SDat , list(nondups[1:25]))
#' 
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedPeptides}}
#' @export
#' @param data - msExperiment
#' @param dups - list as returned by function \code{\link{getListOfMatches}}
analyzeDuplicated = function(data, dups){
  res = NULL
  count = 1
  nameshead = c("p1","p2","cor","rt1","rt2","medianRTDiff","madDiffRT")
  # determine size of output matrix
  tmp = lapply(dups,function(x){d=length(x);return((d*d - d)/2) })
  nrrow = sum( unlist(tmp) )
  cat("nr rows", nrrow , "to process")
  res = matrix( NA , nrow = nrrow ,ncol=length( nameshead ) )
  for(dup in dups)
  {
    duplicated <- dup
    ld =length(duplicated)
    #    print(dup)
    cat("nrdup ", (ld),"\n")
    for(i in 1:ld){
      for(j in i:ld){
        if(i!=j){
          cat("count:",count, " i:",i," j:",j,"\n")
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
  res2= as.data.frame(res[1:(count-1),],stringsAsFactors=FALSE)
  colnames(res2) = nameshead
  res2$cor=as.numeric(res2$cor)
  res2$rt1=as.numeric(res2$rt1)
  res2$rt2=as.numeric(res2$rt2)
  res2$madDiffRT=as.numeric(res2$madDiffRT)
  return(res2)
}

