getwd()
library(data.table)
tmp  =fread("playground/p987_TiMet_6ZTs_sigLib_frgAllReport.txt")
tmp = fread("/Volumes/t1/data/20150730_173210_20150603_TiMet_noDecoys_Report_Frgwise.tsv")

colnames(tmp)
list(
)



colnames(tmp)

unlist(tmp[,1,with=FALSE])

dim(tmp)
for(i in 1: ncol(tmp)){
  x < - (unique(unlist(tmp[,i,with=FALSE])))
  lx <- length(x)
  cat(i , " " , lx , " ", lx/nrow(tmp), "\n")
}

colnames(tmp)
unique(unlist(tmp[,1,with=FALSE]))
