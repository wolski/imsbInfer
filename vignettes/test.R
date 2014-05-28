rm(list=ls())
library(imsbInfer)
SpecLib = ("/home//witold/Analysis/EBhardt//data/E1404301658-sample-SpecLib/feature_alignment_requant.tsv")
xx = fread(SpecLib)

# long running function
data =  transitions2wide(xx)
colnames(data)

##### selecting top 2-n fragments ####
# long running
toptrans = selectTopFragmentsPerPeptide(data,nrt=3)
agrpeptide = aggregatepeptide(toptrans)

# this will read in also the full annotation (which peptide belongs to which protein)
msexp = read2msExperiment(xx)

msexp$Intensity = agrpeptide[,2:dim(agrpeptide)[2]]
rownames(msexp$Intensity) = agrpeptide$transition_group_id


# select top peptides
toppep = selectTopPeptidesPerProtein(msexp ,peptop=3)
boxplot(asinh(toppep$Intensity))

length(toppep$pepinfo$transition_group_id)
# get the transitions belonging to the top peptides

toptrans = toptrans[toppep$pepinfo$transition_group_id]
# make sure that the selection is OK.
xx = aggregatepeptide(toptrans)
toppep$Intensity[xx$transition_group_id[1],1:9] == xx[1,2:10]


