#' get a list of random matches
#' @export
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$PeptideSequence
#' res = getUniquePeptides(rownames)
#' res[1:10]
#' sum(duplicated(res))
#' pairs = randomPeptidePairs(100,res)
#' res2 = analyzeDuplicated( SDat , pairs, maxpep = 100)
#' dim(res2)
#' 
#' @seealso \code{\link{analyzeDuplicated}} for contex
getUniquePeptides = function( rownames )
{
  dups = unique(rownames[which(duplicated(rownames))])
  res = setdiff(rownames,dups)
  res = which( rownames %in% res)
  return(res)
}
#' randomly match peptides
#' @export
#' @examples
#' data(SDat)
#' res = randomPeptidePairs(100,1:dim(SDat)[1])
#' length(res)
#' res[[1]]
#' res2 = analyzeDuplicated( SDat , res)
#' dim(res2)
#' @param number of pairs
#' @seealso \code{\link{analyzeDuplicated}}
randomPeptidePairs = function(n, pepidx){
  res = list(n)
  for(i in 1:n){
    res[[i]] <- sample(pepidx,2)
  }
  return(res)
}
#' get a list of duplicates
#' 
#' @export
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$PeptideSequence
#' res = getListOfMatches(rownames)
#' rownames[res[[1]]]
#' rownames = as.character(SDat$pepinfo$ProteinName)
#' res = getListOfMatches(rownames)
#' lapply(res,length)
#' rownames[res[[1]]]
#' @seealso \code{\link{analyzeDuplicated}} for contex
getListOfMatches = function(  rownames ){
  dups = unique(rownames[which(duplicated(rownames))])
  res = list(length(dups))
  for(i in 1:length(dups)){
    res[[i]]= which(rownames %in% dups[i])
  }
  return(res)
}
#' get a list of duplicates ordered by intensity
#' 
#' @export
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$ProteinName
#' intens = apply(SDat$Intensity,1,mean)
#' res = getListOfMatchesOrderedByIntensity(rownames, apply(SDat$Intensity,1,mean))
#' intens[res[[1]]]
getListOfMatchesOrderedByIntensity = function(  rownames, intensities ){
  dups = unique(rownames[which(duplicated(rownames))])
  res = list(length(dups))
  for(i in 1:length(dups)){
    tmp = which(rownames %in% dups[i])
    xx = intensities[tmp]
    ord = order(xx,decreasing = TRUE)
    res[[i]] = tmp[ord]
  }
  return(res)
}
#' analyse duplicated peptides - same peptide different charge
#'
#' @export
#' @examples
#' data(SDat)
#' res = analyzeDuplicatedPeptides(SDat,countmax= 20)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' #same for proteins
#' rownames = as.character(SDat$pepinfo$ProteinName)
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedProteins}}
#' @param data msExperiment
#' @param countmax maximum number of duplicate comparisons
analyzeDuplicatedPeptides = function(data,countmax = 1000){
  rownames = data$pepinfo$PeptideSequence
  dups=getListOfMatches( rownames )
  analyzeDuplicated( data ,dups, countmax = countmax )
}
#' analyse duplicated protein - same protein id different peptide
#'
#' @export
#' @examples
#' data(SDat)
#' res = analyzeDuplicatedProteins(SDat)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' #same for proteins
#' rownames = as.character(SDat$pepinfo$ProteinName)
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedPeptides}}
#' @seealso  \code{\link{analyzeDuplicated}}
#' @param data an object of class msexperiment
#' @param countmax maximum number of duplicate comparisons
analyzeDuplicatedProteins = function(data,maxpep=3,countmax = 1000){
  rownames = data$pepinfo$ProteinName
  dups=getListOfMatches( rownames )
  analyzeDuplicated( data ,dups , maxpep= maxpep,countmax=countmax)
}
#' analyse duplicated protein - same protein id different peptide take top peptides
#'
#' @export
#' @examples
#' data(SDat)
#' res = analyzeDuplicatedProteins(SDat)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' #same for proteins
#' rownames = as.character(SDat$pepinfo$ProteinName)
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedPeptides}}
#' @seealso  \code{\link{analyzeDuplicated}}
#' @param data an object of class msexperiment
#' @param countmax maximum number of duplicate comparisons
analyzeDuplicatedProteinsTOP = function(data,maxpep=3,countmax = 1000){
  rownames = data$pepinfo$ProteinName
  dups=getListOfMatchesOrderedByIntensity( rownames, apply(data$Intensity,1,mean) )
  analyzeDuplicated( data ,dups , maxpep= maxpep,countmax=countmax)
}

#' analyse duplicated peptides/proteins/lines
#' @examples
#' data(SDat)
#' rownames = SDat$pepinfo$PeptideSequence
#' dups = getListOfMatches(rownames)
#' res = analyzeDuplicated(SDat , dups[1:25])
#' dim(res)
#' res[1,]
#' hist(res$cor)
#' plot(res$medianRTDiff,res$cor)
#' 
#' rownames = SDat$pepinfo$PeptideSequence
#' nondups = getUniquePeptides(rownames)
#' length(nondups)
#' res = analyzeDuplicated(SDat , list(nondups[1:25]))
#' 
#' @seealso \code{\link{getListOfMatches}} \code{\link{analyzeDuplicatedPeptides}}
#' @export
#' @param data - msExperiment
#' @param dups - list as returned by function \code{\link{getListOfMatches}}
#' @param maxpep - limit the number of duplicated peptides
#' @param countmax - limit the total number of comparisons
analyzeDuplicated = function(data, dups, maxpep=3, countmax = 1000, method="pearson"){
  res = NULL
  count = 1
  nameshead = c("p1","p2","cor","rt1","rt2","medianRTDiff","madDiffRT")
  # determine size of output matrix
  tmp = lapply(dups,function(x){d=min(maxpep,length(x));return((d*d - d)/2) })
  nrrow = sum( unlist(tmp) )
  cat("nr rows", nrrow , "to process\n")
  res = matrix( NA , nrow = nrrow ,ncol=length( nameshead ) )
  print(dim(res))
  
  for(dup in dups){
    duplicated <- dup
    ld = min(maxpep, length(duplicated))
    #cat("nrdup ", (ld), " count " , count, "\n")
    if(count > countmax){
      break
    }
    medianrt = apply(data$rt,1,median)
    
    for(i in 1:ld){
      for(j in i:ld){
        if(i!=j){
          cat("count:",count, " i:",i," j:",j,"\n")
          tmp=cor(t(data$Intensity[duplicated[c(i,j)],]), use="pairwise.complete.obs", method=method )
          cors = tmp[upper.tri(tmp)]
          x<-which(upper.tri(tmp),arr.ind=T)
          nams <- rownames(tmp)[x]
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
#' compute various correlation ...
#' @export
compPeptideCorrelations = function(specLib,countmax= 500){
  xxPeptide = analyzeDuplicatedPeptides( specLib , countmax=countmax )
  xxProt = analyzeDuplicatedProteins( specLib , countmax=countmax )
  xxProt2 = analyzeDuplicatedProteinsTOP( specLib , countmax=countmax )
  # look at random peptide pairs
  res = randomPeptidePairs(countmax,1:dim(specLib)[1])
  resRandom = analyzeDuplicated( specLib , res)
  # look at unique peptides
  rownames = specLib$pepinfo$PeptideSequence
  res = getUniquePeptides(rownames)
  pairs = randomPeptidePairs(countmax,res)
  res2 = analyzeDuplicated( specLib , pairs)
  pp = list(peptide = xxPeptide, protein = xxProt, proteinTOP = xxProt2, random = resRandom, unique= res2)
  return(pp)
}
#' plot the output of compPeptideCorrelations
#' @export
plotPairCors=function(res,main="",xlim=c(-1,1))
{
  Dpep = density(res$peptide$cor,na.rm = TRUE)
  Drandom = density(res$random$cor,na.rm = TRUE)
  Dprotein = density(res$protein$cor,na.rm = TRUE)
  DproteinTOP = density(res$proteinTOP$cor,na.rm = TRUE)
  Dunique = density(res$unique$cor,na.rm = TRUE)
  
  maxY = max(max(Dpep$y),max(Drandom$y),max(Dprotein$y),max(DproteinTOP$y),max(Dunique$y))
  
  plot( Dpep, ylim=c(0 , maxY) ,xlim=xlim ,main=main , col=5)
  lines( Drandom, col=2)
  lines( Dprotein, col=3 )
  lines( DproteinTOP, col=6, lty=2)
  lines( Dunique , col=1)

  legend("topleft",legend=c("unique","random","protein","protein Top","peptide"),lty=c(1,1,1,2,1), lwd = c(2,2,2,2,2) , col = c( 1, 2, 3, 6, 5))
  abline(v=0)
}
#' summaryCors
#' @export
summaryCors = function(res){
  xx = cbind(summary(res[[1]]$cor),
         summary(res[[2]]$cor)
         ,summary(res[[3]]$cor)
         ,summary(res[[4]]$cor)
         ,summary(res[[5]]$cor))
  colnames(xx) = names(res)
  return(xx)
}
