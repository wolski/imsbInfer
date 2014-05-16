if(getRversion() >= "3.1.0") utils::suppressForeignCheck("localvariable")
#' order by RT
#' @param obj obj
#' @export
orderByRT = function(obj){
  UseMethod('orderByRT')
}
#' orders rows by retention time of peptide
#' @param experiment obj
#' @export
#@S3method orderByRT msexperiment
#' @examples
#' data(SDat)
#' head(SDat$rt)
#' SDat=orderByRT(SDat)
#' head(SDat$rt)
#' SDat$RT[1:10] 
orderByRT.msexperiment = function(experiment){
  RT = apply(experiment$rt , 1 , median, na.rm=TRUE )
  #order relevant data by retention time
  rto = order(RT)
  experiment$intensity = experiment$intensity[rto,]
  experiment$score = experiment$score[rto,]
  experiment$rt = experiment$rt[rto,]
  experiment$pepinfo = experiment$pepinfo[rto,]
  experiment$RT = RT[rto]
  return(experiment)
}
#' remove decoys
#' @param obj object
#' @export
removeDecoys = function(obj,...){
  UseMethod('removeDecoys')
}
#' remove decoys from msexperiment
#' @param data msexperiment
#' @export
#@S3method removeDecoys msexperiment
#' 
#' @examples
#' data(feature_alignment_requant)
#' SDat = read2msExperiment(feature_alignment_requant)
#' dim(SDat)
#' SDat = removeDecoys(SDat)
#' dim(SDat)
removeDecoys.msexperiment = function(data){
  return(subset(data,!data$pepinfo$decoy))
}
#' removes unwanted RT ranges
#' @param obj object
#' @export
keepRTRange <- function(obj, ...){
  UseMethod('keepRTRange')
}
#' filter data given RT range
#'
#' @param data msexperiment
#' @param rtrange rt range to keep default 1000 - 7000
#' @export
#@S3method keepRTRange msexperiment
#' @examples
#' data(feature_alignment_requant)
#' SDat = read2msExperiment(feature_alignment_requant)
#' dim(SDat)
#' SD2 = keepRTRange(SDat,c(3000,4000))
#' dim(SD2)
keepRTRange.msexperiment<- function(data,rtrange=c(1000,7000)){
  RT = apply( data$rt , 1 , median )
  idx = RT > rtrange[1] &  RT < rtrange[2]
  return( subset( data , idx ) )
}
#' @export
subset <- function(obj, ...){
  UseMethod('subset')
}
#' subset data given idx vector
#' @aliases subset
#' @param data msexperiment
#' @param idx row indices to keep or of TRUE, FALSE vector
#' @export
#' @examples
#' data(SDat)
#' # keep only peakgroups with mass less than 800
#' mass = apply(SDat$mz,1,median)
#' dim(SDat)
#' SDatr = subset(SDat, mass < 800)
#' dim(SDatr)
subset.msexperiment<- function(data,idx){
  data$intensity = data$intensity[idx,]
  data$score = data$score[idx,]
  data$rt = data$rt[idx,]
  data$pepinfo = data$pepinfo[idx,]
  if(!is.null(data$RT)){
    data$RT = data$RT[idx]
  }
  return(data)
}
#' convert LF 2 wide format
#' @param aligtable output of RT alignment tool (long format)
#' @export
#' @examples
#' data(feature_alignment_requant) 
#' SDat =  convertLF2Wideformat(  feature_alignment_requant  ) 
#'
convertLF2Wideformat=function(aligtable){
  
  aligtable = as.data.table(aligtable)
  idx = grep("m*z" , colnames(aligtable) )
  setnames(aligtable,  idx, "mz")
  #fix orignal filename
  aligtable$align_origfilename = basename(as.character(aligtable$align_origfilename))
  protmapping = aligtable[,list(transition_group_id,ProteinName)]
  setkey(protmapping,transition_group_id)
  
  #fix missing keys by adding NA
  setkey(aligtable,transition_group_id,align_origfilename)
  aligtable = aligtable[CJ(unique(transition_group_id), unique(align_origfilename))]
  # extract intensitiies
  ints = aligtable[,as.list(Intensity),by=transition_group_id]
  setnames(ints,c("transition_group_id",paste("Intensity", unique(aligtable$align_origfilename ), sep="_")))
  #extract rt's
  rt = aligtable[,as.list(RT),by=transition_group_id]
  setnames(rt,c("Peptide",paste("RT", unique(aligtable$align_origfilename ), sep="_")))
  #extract scor's
  score = aligtable[,as.list(m_score),by=transition_group_id]
  setnames(score,c("Peptide",paste("score", unique(aligtable$align_origfilename ), sep="_")))
  #extract masses
  mz = aligtable[,as.list(mz),by=transition_group_id]
  setnames(mz,c("Peptide",paste("mz", unique(aligtable$align_origfilename ), sep="_")))
  
  return(list(protmapping=unique(protmapping), ints=ints, rt=rt,score=score,mz=mz, aligtable=aligtable))
}
#' converts to msExperiment (exported more for debugging purpose)
#' @param data output of convertLF2Wideformat
#' @export
#' @examples
#' data(feature_alignment_requant) 
#' SDat=convert2msExperiment(convertLF2Wideformat(feature_alignment_requant))
#' 
convert2msExperiment = function(data){
  nams=colnames(data$ints)
  # prepare the column names
  nams = sub("Intensity_","",nams)
  nams = sub("_with_dscore.*","",nams)
  nams = sub("_all_peakgroups.*","",nams)
  setnames(data$ints,nams)
  setnames(data$score,nams)
  setnames(data$rt, nams)
  setnames(data$mz, nams)
  
  #setnames(data$protmapping,nams)
  nams = data$protmapping$transition_group_id
  idxDecoy<-grep("^DECOY\\_",nams)
  decoy = rep(FALSE,length(nams))
  decoy[idxDecoy] <- TRUE
  
  nams = sub("^DECOY\\_","",nams)
  
  pepinfo=split2table(nams,split="\\_|\\/")
  pepinfo = pepinfo[,-dim(pepinfo)[2]]
  colnames(pepinfo)=c("id","sequence","charge")
  pepinfo = as.data.frame(pepinfo,stringsAsFactors = FALSE)
  
  nametable = data.frame(transition_group_id=as.character(data$protmapping$transition_group_id),
                         pepinfo,
                         decoy=decoy,
                         seqchid = paste(pepinfo$sequence, pepinfo$charge, sep="_"),
                         ProteinName = as.character(data$protmapping$ProteinName),
                         stringsAsFactors = FALSE
                         )

  #nametable  = nametable[order(nametable$transition_group_id),]
  nrcol = dim(data$ints)[2]
  SwathDat = list(intensity = as.matrix(data$ints[,2:nrcol,with=F] ) ,
                  score = as.matrix(data$score[,2:nrcol,with=F] ) ,
                  rt  = as.matrix(data$rt[,2:nrcol,with=F]) ,
                  mz  = as.matrix(data$mz[,2:nrcol,with=F]) ,
                  pepinfo = nametable
                  )
  
  rownames(SwathDat$intensity) = data$ints$transition_group_id
  rownames(SwathDat$score) = data$ints$transition_group_id
  rownames(SwathDat$rt) = data$ints$transition_group_id
  rownames(SwathDat$mz) = data$ints$transition_group_id
  
  colnames(SwathDat$intensity) = colnames(data$ints)[2:nrcol]
  colnames(SwathDat$score) = colnames(data$ints)[2:nrcol]
  colnames(SwathDat$rt) = colnames(data$ints)[2:nrcol]
  colnames(SwathDat$mz) = colnames(data$ints)[2:nrcol]
  
  class(SwathDat) = "msexperiment"
  return(SwathDat)
}

#' @export
read2msExperiment=function(obj,...){
  UseMethod('read2msExperiment')
}
#' read 2 feature alginer long format and generate an msexperiment
#' @aliases read2msExperiment
#' @param filename filename to read (can't be compressed)
#' @return msexperiment
#' @export
read2msExperiment.default=function(filename){
  print("read2msExperiment.default")
  aligtable=fread(filename)
  data = convertLF2Wideformat(aligtable)
  data = convert2msExperiment(data)
  return(data)
}
#' convert data.frame 2 msexperiment
#' @param data a data.frame in long format
#' @return msexperiment
#' @export
#' @examples
#' data(feature_alignment_requant)
#' SDat = read2msExperiment(feature_alignment_requant)
#' dim(SDat)
read2msExperiment.data.frame=function(data){
  data = convertLF2Wideformat(data)
  data = convert2msExperiment(data)
  return(data)
}
#' dimension
#'
#' @export
#@S3method dim msexperiment
dim.msexperiment<-function(x){
  return(dim(x$intensity))
}
#' show colnames (it does'nt let you set the columns)
#' 
#' @examples
#' data(feature_alignment_requant)
#' SDat = read2msExperiment(feature_alignment_requant)
#' colnames(SDat)
#' @export
#@S3method colnames msexperiment
colnames.msexperiment<-function(x){
  return(colnames(x$intensity))
}
#' allows to set colnames for all the matrices in msexperiment
#'
#' @export
mycolnames  = function(data,newNames)
{
  for(i in 1:(length(data)-1))
  {
    colnames(data[[i]]) = newNames
  }
  return(data)
}
#' reorders all the matrices by the columnnames
#' @param msexperiment to reorder
#' @param ord new ordering - if null then use ordering by column names
#' @export
reordercolumns  = function( data, ord = NULL )
{
  for(i in 1:(length(data)-1))
  {
    if(length(ord)==0){
      ord= order(colnames(data[[i]]))
    }
    data[[i]] = data[[i]][,ord]
  }
  return(data)
}



