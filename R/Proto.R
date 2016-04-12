rm(list=ls())
gc()

SpecLib = "C:/Users/wewol/Google Drive/tissuecomparison/OpenSWATH/BAT_19strains/data/E1603291025_feature_alignment_requant.tsv"

library(readr)
data <- read_tsv(SpecLib,col_names = TRUE)
#tmp <-read.csv(SpecLib, stringsAsFactors = F,sep="\t",header=T)

source("R/openSwath.R")
head(data)

glimpse(data)
data2 <- prepareDF(data)
head(data2)
data3 <- prepareOpenSwathData(data2)
head(data3)
dim(transids)

library(reshape2)
transIntensities = dcast(data3, transition_group_id + ion_type + frag_charge + pep_sequence ~ align_origfilename , value.var="aggr_Peak_Area")

tmp <-split2table(as.character(transIntensities$aggr_Fragment_Annotation) , split = "_")
colnames(tmp)<-c("frag_id", "ion_type","fragCharge","pepSequence","pepCharge")


head(transIntensities)


mscore = dcast(data2, transition_group_id  ~ align_origfilename , value.var="m_score")
rt = dcast(data2, transition_group_id ~ align_origfilename , value.var="RT")
mz = dcast(data2, transition_group_id ~ align_origfilename , value.var="mz")
pepIntensities = dcast(data2,  transition_group_id ~ align_origfilename , value.var="Intensity")

head(data3)
split2table(data)

