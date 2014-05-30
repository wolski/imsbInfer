rm(list=ls())
library(imsbInfer)

SpecLib = ("/home//witold/Analysis/EBhardt//data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
data = fread(SpecLib)

colnames(data)
# keep only non decoy
sum(data$decoy==1)
data = data[data$decoy==0]

# nr of transition per peptide
nrt = 3
# nr of peptides per protein
peptop = 3

# 
msexp = loadTransitonsMSExperiment(data,nrt=nrt,peptop=peptop)
table(table(msexp$pepinfo$ProteinName))


msexpS = read2msExperiment(data)
length(unique(msexpS$pepinfo$ProteinName))

xx = selectTopPeptidesPerProtein(msexpS,peptop=3)
length(unique(xx$pepinfo$ProteinName))






tmp = convert2MSstats(msrob)

head(tmp)

msexp2 = aggregatePeptides(msexp)
msexp3 = aggregateProtein(msexp)


#####
# reading tiannan dataset : 6_049_920 rows
library(imsbInfer)
SpecLib = ("/home//witold/Analysis/tiannan/data/feature_alignment_requant.tsv")
data = fread(SpecLib)
data = data[data$decoy==0]
msexp = loadTransitonsMSExperiment(data,nrt=3,peptop=3)


