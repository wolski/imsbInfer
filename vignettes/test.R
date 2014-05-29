rm(list=ls())
library(imsbInfer)


SpecLib = ("/home//witold/Analysis/EBhardt//data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
data = fread(SpecLib)

nrt = 3
peptop = 3

msexp = loadTransitonsMSExperiment(data,nrt=nrt,peptop=peptop)
boxplot(asinh(msexp$Intensity))


SpecLib = ("/home//witold/Analysis/LightSwath/MSstatsanalysis/feature_alignment_requant.tsv")
data = fread(SpecLib)

nrt = 3
peptop = 3

obj = data
tmp = basename(obj$align_origfilename)
tmp[1:10]
tmp = sub("_SW_with_dscore.csv","",tmp)

obj$align_origfilename = paste(tmp , obj$align_runid , sep="_")
obj$align_origfilename[1:10]
data =  transitions2wide(obj)
toptrans = selectTopFragmentsPerPeptide(data,nrt=3)

dim(data)
colnames(data)
msexp = loadTransitonsMSExperiment(data,nrt=nrt,peptop=peptop)
boxplot(asinh(msexp$Intensity))

