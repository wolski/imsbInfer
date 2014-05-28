#' extract transition intensities from input file. WARNING - slow running function
#' @description
#' the input table required fields are:
#' transition_group_id, 
#' aggr_Peak_Area,
#' aggr_Fragment_Annotion,
#' align_origfilename, 
#' @export
#' @examples
#' data(feature_alignment_requant)
#' transitions2wide(feature_alignment_requant)
transitions2wide = function(far){
  #far = feature_alignment_requant
  ids = as.character(far$transition_group_id)
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  orig = basename(as.character(far$align_origfilename))
  orig = sub("_SW_with_dscore.csv","",orig)
  
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";"))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";"))})
  
  # prepare output
  lx = length(transids)
  idss = vector(length=lx,mode="list")
  origf = vector(length=lx,mode="list")
  
  # extend
  for(i in 1:lx){
    l =  length(transids[[i]])
    idss[[i]] = rep(ids[i],l) # transition group ids
    origf[[i]] = rep(orig[i],l) # orig file name
  }
  library(reshape2)
  # extend to wide format
  tmp = data.frame(transition_group_id = as.character(unlist(idss)),
                   align_origfilename = as.character(unlist(origf)),
                   aggr_Peak_Area = as.numeric(unlist(transints)),
                   aggr_Fragment_Annotation =  as.character(unlist( transids) ) )
  #got wide format here...
  data = dcast(tmp, transition_group_id + aggr_Fragment_Annotation ~ align_origfilename , value.var="aggr_Peak_Area")
  return(data)
}
#' selecting top transtions WARNING - slow running function
#' @export
#' @param data - data.frame or table with at least 2 columns:  transition_group_id, aggr_Fragment_Annotation (unique keys)
#' @param nrt - top transitions to select
#' @examples
#' data(feature_alignment_requant)
#' tmp = transitions2wide(feature_alignment_requant)
#' xx = selectTopFragmentsPerPeptide(tmp)
selectTopFragmentsPerPeptide = function(data, nrt = 2  ){
  #compute median and create table with id's
  medxx = apply(data[,3:dim(data)[2]],1,median)
  xxmex= cbind( data[,c("transition_group_id","aggr_Fragment_Annotation")] , medxx)
  tmpdt = data.table(xxmex)
  # fixing column types
  tmpdt <- tmpdt[, transition_group_id:=as.character(transition_group_id)]
  tmpdt <- tmpdt[, aggr_Fragment_Annotation:=as.character(aggr_Fragment_Annotation)]
  
  setkey(tmpdt , transition_group_id,aggr_Fragment_Annotation)
  
  # we are going to select the top peptides for each transition group.
  transgroupid = unique( as.character(tmpdt$transition_group_id) )
  
  # prepare output matrix
  res <- matrix("",nrow=length(transgroupid)*nrt,ncol=3)
  dim(res)
  start = 1
  end = 0
  for(i in 1:length(transgroupid)){
    tmp<-tmpdt[transgroupid[i],]
    xx=which(order(as.numeric(tmp$medxx),decreasing=TRUE) < (nrt+1))
    end = start + (length(xx)-1)
    xbla = as.matrix(tmp[xx,])
    res[start:end,] = xbla
    start = end + 1
  }
  res = res[1:end,]
  colnames(res) = colnames(tmpdt)
  res = data.table(res)
  data = data.table(data)
  setkey(res,transition_group_id,aggr_Fragment_Annotation)
  setkey(data,transition_group_id,aggr_Fragment_Annotation)
  res = data[res]
  # drop last row which is the median
  res = res[,-dim(res)[2],with=FALSE]
  return(res)
}

#' aggregate peptide from transtion
#' @export
#' @examples
#' NULL
aggregatepeptide=function(toptrans, func = sum){
  toptransvals=toptrans[,3:dim(toptrans)[2],with=F]
  toptransvals = lapply(toptransvals,as.numeric)
  agregatepeptide = aggregate(toptransvals,by=list(toptrans$transition_group_id),func)
  colnames(agregatepeptide)[1] = "transition_group_id"
  # select data values
  agregatepeptide$transition_group_id = as.character(agregatepeptide$transition_group_id)
  return(agregatepeptide)
}

#' this function selects the top x peptides / protein
#' @param newprot data.frame with 
#' @export
selectTopPeptidesPerProtein = function(msexp, peptop = 3){
  #newprot = merge(msexp$pepinfo[,c("transition_group_id","ProteinName")],agrpeptide,by.x="transition_group_id",by.y="transition_group_id")
  
  #compute median and create table with id's
  medxx = apply(msexp$Intensity , 1,median)
  xxmex = cbind( msexp$pepinfo[,c("transition_group_id","ProteinName")] , medxx)
  
  tmpdt = data.table(xxmex)
  ## fixing column types
  tmpdt <- tmpdt[, transition_group_id:=as.character(transition_group_id)]
  tmpdt <- tmpdt[, ProteinName:=as.character(ProteinName)]
  setkey(tmpdt , ProteinName,transition_group_id)
  
  # we are going to select the top peptides for each protein.
  proteinname = unique( as.character(tmpdt$ProteinName) )
  # prepare output matrix
  lx = length(proteinname)
  res <- matrix("",nrow=lx*peptop,ncol=3)
  start = 1
  end = 0
  dim(res)
  dim(tmpdt)
  for(i in 1:lx){
    tmp<-tmpdt[proteinname[i],]
    xx=which(order(as.numeric(tmp$medxx),decreasing=TRUE) < (peptop+1))
    end = start + (length(xx)-1)
    xbla = as.matrix(tmp[xx,])
    res[start:end,] = xbla
    start = end + 1
  }
  
  res = res[1:end,]
  colnames(res) = c("ProteinName","transition_group_id","medxx")

  res = data.table(res)
  res = subset(msexp, res$transition_group_id)
  return(res)
}
