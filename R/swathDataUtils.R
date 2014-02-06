#' filter data given RT range
#'
#' @export
removeRTRanges<- function(data,rtrange=c(1000,7000)){
  idx = SwathDat$RT > rtrange[1] &  SwathDat$RT < rtrange[2]
  SwathDat$peplev = SwathDat$peplev[idx,]
  SwathDat$score = SwathDat$score[idx,]
  SwathDat$newid2 = SwathDat$newid2[idx ]
  SwathDat$newid = SwathDat$newid[idx ]
  SwathDat$protnam = SwathDat$protnam[idx]
  SwathDat$RT = SwathDat$RT[idx ]
}
#' read 2 matrix export return msExperiment
#' 
#' @export
read2matrixExport=function(filename){
  swathpep = read.table(filename,header=T,stringsAsFactors=FALSE)
  #split into intensities and scores
  swathpeplev = swathpep[,grep("Intensity",colnames(swathpep))]
  swathscore = swathpep[,grep("^score",colnames(swathpep))]
  
  # prepare the column names
  nams = sub("Intensity_","",colnames(swathpeplev))
  nams = sub("_with_dscore.*","",nams)
  nams = sub("_all_peakgroups.*","",nams)
  
  colnames(swathscore) = nams
  colnames(swathpeplev) = nams

  #remove decoys
  decoys = grep("DECOY",swathpep[,1])
  swathpeplev = swathpeplev[-decoys,]
  protnam = swathpep[-decoys,1:2]
  swathscore = swathscore[ -decoys, ]
  swathpepRT = swathpep$RT_mean[ -decoys ]
  
  # prepare peptide information
  pepseq = protnam[,1]
  pepinfo=split2table(pepseq,split="\\_")
  
  if(dim(swathpeplev)[1] != length(swathpepRT) || 
  dim(swathpeplev)[1] != dim(protnam)[1] ||
  dim(swathpeplev)[1] != dim(pepinfo) ||
  dim(swathpeplev)[1] != dim(swathscore)[1]){
    stop("dimension error")
  }
  
  res = list(filename = filename,
                  pepinfo=pepinfo,
                  protnam = protnam[,2],
                  RT = swathpepRT,
                  peakscore = swathscore,
                  peplev = swathpeplev 
  )
  class(res) <- "msexperiment"
  return(res)
}