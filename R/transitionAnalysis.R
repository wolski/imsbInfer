#' checks input file for required columns - removes nonrequired columns from data.frame
#' @param df ouptut fo OpenSwath feature aligner
#' @export
#' @examples
#' library(imsbInfer)
#' data(feature_alignment_requant)
#' colnames(feature_alignment_requant)
#' df = feature_alignment_requant
#' df = prepareDF(feature_alignment_requant)
#' colnames(df)
#' length(unique(df$ProteinName))
prepareDF <- function(df){
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  required = c("transition_group_id","align_origfilename","RT",
               "mz","Intensity","ProteinName","m_score","aggr_Fragment_Annotation","aggr_Peak_Area")
  df  = data.table(df)
  
  setnames(df, sub("m.z","mz",colnames(df)))
  x=match(required,colnames(df))
  stopifnot(required == colnames(df)[x])
  df = df[,x,with=FALSE]
  df <- df[order(df$transition_group_id),]
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(df)
}
#' extract transition intensities from input file. WARNING - slow running function
#' @description
#' the input table required fields are:
#' transition_group_id, 
#' aggr_Peak_Area,
#' aggr_Fragment_Annotion,
#' align_origfilename,
#' 
#' @param far convert to wide format 
#' @export
#' @examples
#' data(feature_alignment_requant)
#' df = prepareDF(feature_alignment_requant)
#' tmp = transitions2wide(df)
#' names(df)
#' dim(feature_alignment_requant)
#' dim(tmp)[1]/dim(feature_alignment_requant)[1]*3
transitions2wide <- function(far){
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  #far = feature_alignment_requant
  ids = as.character(far$transition_group_id)
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  orig = basename(as.character(far$align_origfilename))
  orig = sub("_with_dscore.*","",orig)
  
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";"))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";"))})
    # prepare output
  lx = length(transids)
  transids <<-transids
  idss = vector(length=lx,mode="list")
  origf = vector(length=lx,mode="list")
  
  # extend
  for(i in 1:lx){
    l =  length(transids[[i]])
    idss[[i]] = rep(ids[i],l) # transition group ids
    origf[[i]] = rep(orig[i],l) # orig file name
  }
  # extend to wide format
  tmp = data.frame(transition_group_id = as.character(unlist(idss)),
                   align_origfilename = as.character(unlist(origf)),
                   aggr_Peak_Area = as.numeric(unlist(transints)),
                   aggr_Fragment_Annotation =  as.character(unlist( transids) ) )
  #got wide format here...
    data = dcast(tmp, transition_group_id + aggr_Fragment_Annotation ~ align_origfilename , value.var="aggr_Peak_Area")
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(data)
}
#' selecting top transtions WARNING - slow running function
#' @export
#' @param data - data.frame or table with at least 2 columns:  transition_group_id, aggr_Fragment_Annotation (unique keys)
#' @param nrt - top transitions to select
#' @examples
#' data(feature_alignment_requant)
#' table(table(feature_alignment_requant$transition_group_id))
#' tmp = transitions2wide(feature_alignment_requant)
#' # how many peptides with how many transitions
#' x1 = table(table(tmp$transition_group_id))
#' xx = selectTopFragmentsPerPeptide(tmp,nrt=7)
#' x2 = table(table(xx$transition_group_id))
#' # because with nrt 7 you do not remove any transitions
#' stopifnot(x1 == x2)
#' xx = selectTopFragmentsPerPeptide(tmp,nrt=2)
#' # all peptides must have 2 transitions
#' stopifnot(names(table(table(xx$transition_group_id))) == "2")
selectTopFragmentsPerPeptide = function(data, nrt = 2  ){
  loccoll = Sys.getlocale("LC_COLLATE")
  
  Sys.setlocale("LC_COLLATE", "C")
  
  #compute median and create table with id's
  medxx = apply(data[,3:ncol(data)],1, mean ,na.rm=TRUE)
  xxmex= cbind( data[,c("transition_group_id","aggr_Fragment_Annotation")] , medxx)
  tmpdt = data.table(xxmex)
  # fixing column types
  tmpdt <- tmpdt[, transition_group_id:=as.character(transition_group_id)]
  tmpdt <- tmpdt[, aggr_Fragment_Annotation:=as.character(aggr_Fragment_Annotation)]
  
  setkey(tmpdt , transition_group_id,aggr_Fragment_Annotation)
  
  # we are going to select the top transitions for each transition group.
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
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(res)
}
#' aggregate peptide from transtion
#' @export
#' @examples
#' data(feature_alignment_requant)
#' tmp = transitions2wide(feature_alignment_requant)
#' xx = selectTopFragmentsPerPeptide(tmp)
#' aggr = .aggregatepeptide(xx)
#' dim(aggr)
.aggregatepeptide <- function(toptrans, func = sum){
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  # this is a compatibility hack since data.table does not work on windows...
  toptrans = as.data.frame(toptrans)
  toptransvals=toptrans[,3:dim(toptrans)[2]]
  toptransvals = lapply(toptransvals,as.numeric)
  agregatepeptide = aggregate(toptransvals,by=list(toptrans$transition_group_id),func)
  colnames(agregatepeptide)[1] = "transition_group_id"
  # select data values
  agregatepeptide$transition_group_id = as.character(agregatepeptide$transition_group_id)
  Sys.setlocale("LC_COLLATE", loccoll)
  
  return(agregatepeptide)
}

#' this function selects the top x peptides / protein
#' @param msexp data.frame with
#' @param peptop how many top peptides
#' @export
#' @examples
#' data(SDat)
#' rownames(SDat$pepinfo)
#' table(table(SDat$pepinfo$ProteinName))
#' x = selectTopPeptidesPerProtein(SDat,peptop=3)
#' table(table(x$pepinfo$ProteinName))
#' stopifnot(rownames(x$pepinfo)[1:10]==rownames(x$Intensity)[1:10])
#' stopifnot( length(unique(SDat$pepinfo$ProteinName)) == length(unique(x$pepinfo$ProteinName)) )
selectTopPeptidesPerProtein <- function(msexp, peptop = 3){
  loccoll = Sys.getlocale("LC_COLLATE")
  Sys.setlocale("LC_COLLATE", "C")
  
  #newprot = merge(msexp$pepinfo[,c("transition_group_id","ProteinName")],agrpeptide,by.x="transition_group_id",by.y="transition_group_id")
  #compute median and create table with id's
  medxx = apply(msexp$Intensity , 1,mean,na.rm=TRUE)
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
  for(i in 1:lx){
    tmp<-tmpdt[proteinname[i],]
    xx=which(order(as.numeric(tmp$medxx),decreasing=TRUE) < (peptop+1))
    end = start + (length(xx)-1)
    xbla = as.matrix(tmp[xx,])
    res[start:end,] = xbla
    start = end + 1
  }
  
  res = res[1:end,]
  head(res)
  colnames(res) = c("transition_group_id","ProteinName","medxx")
  res = data.table(res)
  dim(res)
  length(unique(msexp$pepinfo$ProteinName))
  length(unique(res$ProteinName))
  
  res2 = subset(msexp, match(res$transition_group_id, rownames(msexp$Intensity)))
  head(msexp$pepinfo$transition_group_id)
  rownames(msexp$Intensity)
  
  head(res2$pepinfo)
  
  #postconditions
  stopifnot(res2$pepinfo$ProteinName==res$ProteinName)
  stopifnot(length(unique(msexp$pepinfo$ProteinName))==length(unique(res$ProteinName)))
  
  Sys.setlocale("LC_COLLATE", loccoll)
  
  
  return(res2)
}

