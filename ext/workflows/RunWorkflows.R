library(compiler)
enableJIT(2)

rm(list=ls())
SpecLib = ("data/E1503161309_feature_alignment_requant.tsv")

Sweave("analyseSWATH.Rnw")
tools::texi2dvi("analyseSWATH.tex",pdf=TRUE)


########set parameters
rm(list=ls())
SpecLib = ("data/E1503161309_feature_alignment_requant.tsv")
Samples2Remove = "ywu_K150219_007_SW"
FileSampleMapping = "data/mappingNew.tsv"

#ouptups
MSSTATS_INPUT="output/MSStatsInputFile.tsv"
MSEXPERIMENT="output/msexporiTrans.rda"

#parameters - map ms stats Conditions to Tissue
Condition="Tissue"
BioReplicate="Strain"

#filter transitions/peptides by mscore
scorethresh=0.005
fraction=0.25

#how many transitions to select
maxNRTransitions = 5
maxNRPeptides = 10

Sweave("PrepareForMSStats.Rnw")
tools::texi2dvi("PrepareForMSStats.tex",pdf=TRUE)


#################

rm(list=ls())
library(imsbInfer)
TOPPeptidesPerProtein = 3
# this file is generated using the PrepareForMSStatsEH2 script...
RAW_OUTPUT="output/ProteinQuant3TOP3TOPRAW.txt"
TRANSFORMED_OUTPUT="output/ProteinQuant3TOP3TOPRAWTransformed.txt"

Sweave("aggregateProtein.Rnw")
tools::texi2dvi("aggregateProtein.tex",pdf=TRUE)

#################
MSSTATS_INPUT="output/MSStatsInputFile.tsv"
MSSTATS_QUANT_OUTPUT="output/MsStatsQuantResults.txt"

Sweave("MSSatsQuantification.Rnw")
tools::texi2dvi("MSSatsQuantification.tex",pdf=TRUE)
