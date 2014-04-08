#xx = read.table("/home/witold/prog/imsbInfer/data/feature_alignment_requant.csv",header=TRUE,sep=",")
#write.table(xx,file="/home/witold/prog/imsbInfer/data/feature_alignment_requant.csv",sep=";")
#' read 2 matrix export return msExperiment (must contain score columns)
#' @param filename output of feature aligner
#' @param decoy.rm
#' @export
deprec.read2msExperiment=function(filename,decoy.rm=TRUE){
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
  
  ord = order(nams)
  swathscore = swathscore[,ord]
  swathpeplev = swathpeplev[,ord]
  
  #swathpeplev = swathpeplev[-decoys,]
  protnam = swathpep[,1:2]
  swathpepRT = swathpep$RT_mean
  
  #remove decoys
  if(decoy.rm){
    decoys = grep("DECOY",swathpep[,1])
    swathpeplev = swathpeplev[-decoys,]
    protnam = swathpep[-decoys,]
    swathscore = swathscore[ -decoys, ]
    swathpepRT = swathpepRT[ -decoys ]
  }
  
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
