rm(list=ls())
library(imsbInfer)


SpecLib = ("/home//witold/Analysis/EBhardt//data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
data = fread(SpecLib)

colnames(data)
sum(data$decoy==1)
data = data[data$decoy==0]
dim(data)
nrt = 3
peptop = 3

msexp = loadTransitonsMSExperiment(data,nrt=nrt,peptop=peptop)
boxplot(asinh(msexp$Intensity))

# do normalization

tmp = convert2MSstats(msexp)

msexp2 = aggregatePeptide(msexp)
dim(msexp2)
head(msexp2$pepinfo)

msexp3 = aggregateProtein(msexp)


#####
# reading tiannan dataset : 6_049_920 rows
library(imsbInfer)
SpecLib = ("/home//witold/Analysis/tiannan/data/feature_alignment_requant.tsv")
data = fread(SpecLib)
data = data[data$decoy==0]
msexp = loadTransitonsMSExperiment(data,nrt=3,peptop=3)


