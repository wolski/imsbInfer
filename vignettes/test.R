library(imsbInfer)
SpecLib = ("/home//witold/Analysis/EBhardt//data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
xx = fread(SpecLib)

colnames(xx)
library(imsbInfer)
data(feature_alignment_requant)
dim(feature_alignment_requant)


far =feature_alignment_requant


transitions2wide = function(far){
  ids = as.character(far$transition_group_id)
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  orig = basename(as.character(far$align_origfilename))
  
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";"))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";"))})
  
  # prapre output
  lx = length(transids)
  idss = list(lx)
  origf = list(lx)
  
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

data = transitions2wide(xx)

##### selecting top 2-n transtions ####
library(data.table)

medxx = apply(data[,3:dim(data)[2]],1,median)
xxmex= cbind( data[1:2] , medxx)
tmpdt = data.table(xxmex)
# fixing column types
tmpdt <- tmpdt[, transition_group_id:=as.character(transition_group_id)]
tmpdt <- tmpdt[, aggr_Fragment_Annotation:=as.character(aggr_Fragment_Annotation)]

setkey(tmpdt , transition_group_id,aggr_Fragment_Annotation)


transid = unique( as.character(tmpdt$transition_group_id) )
length( transid )

## res = vector(length(transid),mode="list")
nrt=3
res <- matrix("",nrow=length(transid)*nrt,ncol=3)
dim(res)
start = 1
end = 0
for(i in 1:length(transid)){
  tmp<-tmpdt[transid[i],]
  xx=which(order(as.numeric(tmp$medxx),decreasing=TRUE) < (nrt+1))
  end = start + (length(xx)-1)
  cat("start",start,"end",end,"\n")
  xbla = as.matrix(tmp[xx,])
  res[start:end,] = xbla
  start = end + 1
}
res = res[1:end,]

colnames(res) = colnames(tmpdt)
dim(res)
head(res)

res = data.table(res)
res


