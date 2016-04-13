rm(list=ls())
gc()


SpecLib = "/data/E1603291025_feature_alignment_requant.tsv"

library(readr)
#tmp <-read.csv(SpecLib, stringsAsFactors = F,sep="\t",header=T)


source("R/openSwath.R")

data2 <- prepareDF(data)
data3 <- prepareOpenSwathData(data2)




dim(precData)



head(transScore)
head(transRT)
head(transMz)
