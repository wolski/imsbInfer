library(imsbInfer)
SpecLib = ("/home//witold/Analysis/EBhardt//data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
xx = fread(SpecLib)

colnames(xx)


library(imsbInfer)
data(feature_alignment_requant)

dim(feature_alignment_requant)
far =feature_alignment_requant
far[13,grep("aggr",names(far))]


ids = as.character(far$transition_group_id)
apa = as.character(far$aggr_Peak_Area)
afa = as.character(far$aggr_Fragment_Annotation)

orig = basename(as.character(far$align_origfilename))

transints = lapply(apa,function(x){unlist(strsplit(x,";"))})
transids = lapply(afa,function(x){unlist(strsplit(x,";"))})

idss = list(length(transids))
origf = list(length(transids))

for(i in 1:length(transids)){
  l=  length(transids[[i]])
  idss[[i]] = rep(ids[i],l) 
  origf[[i]] = rep(orig[i],l)
}

library(reshape2)
library(data.table)

tmp = data.frame(transition_group_id = as.character(unlist(idss)),
                 align_origfilename = as.character(unlist(origf)),
                 aggr_Peak_Area = as.numeric(unlist(transints)),
                 aggr_Fragment_Annotation =  as.character(unlist( transids) ) )
#got wide format here...
data = dcast(tmp, transition_group_id + aggr_Fragment_Annotation ~ align_origfilename , value.var="aggr_Peak_Area")
dim(data)

medxx = apply(xxvalues,1,median)
xxmex= cbind( xx[1:2] , medxx)
tmpdt = data.table(xxmex)
tmpdt <- tmpdt[, transition_group_id:=as.character(transition_group_id)]
tmpdt <- tmpdt[, aggr_Fragment_Annotation:=as.character(aggr_Fragment_Annotation)]
setkey(tmpdt , transition_group_id,aggr_Fragment_Annotation)


transid = unique(as.character(tmpdt$transition_group_id))
length(transid)

#res <- vector(length(transid),mode="list")
nrt=2
res <- matrix("",nrow=length(transid)*nrt,ncol=3)
dim(res)
for(i in 1:length(transid)){
  tmp<-tmpdt[transid[i],]
  xx=which(ord(tmp$medxx) < 3)
  end = nrt*i
  start = (end-nrt+1)
  cat("start",start,"end",end,"\n")
  print(dim(res))
  xbla = as.matrix(tmp[xx,])
  res[start:end,] = xbla
}

colanems(res) = colnames(tmpdt)


xx = unlist(res)
xx[1:12]
xx = matrix(xx,ncol=3,byrow=F)
head(xx)

