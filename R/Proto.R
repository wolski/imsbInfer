rm(list=ls())
gc()

SpecLib = "C:/Users/wewol/Google Drive/tissuecomparison/OpenSWATH/BAT_19strains/data/E1603291025_feature_alignment_requant.tsv"

library(readr)
data <- read_tsv(SpecLib,col_names = TRUE)
#tmp <-read.csv(SpecLib, stringsAsFactors = F,sep="\t",header=T)


source("R/openSwath.R")

data2 <- prepareDF(data)
data3 <- prepareOpenSwathData(data2)




dim(precData)

library(reshape2)

transIntensities = dcast(data4, ModifiedSequence + z + IonType + FragZ + Sequence ~ Filename , value.var="Intensity")

transScore = dcast(precData, ModifiedSequence + z + Sequence ~ Filename , value.var="Score")
transRT = dcast(precData, ModifiedSequence + z + Sequence  ~ Filename , value.var="RT")
transMz = dcast(precData, ModifiedSequence + z + Sequence  ~ Filename , value.var="mz")

head(transScore)
head(transRT)
head(transMz)
