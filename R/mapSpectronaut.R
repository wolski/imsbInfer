.SpectronautPrecursorDefsMapping <- list("FileName"="Run",
                                    "ProteinName"="ProteinName",
                                    "Decoy"="decoy",
                                    "StrippedSequence"="EG.StrippedSequence",
                                    "ModifiedSequence"="PeptideSequence",
                                    #"IsotopeLabelType"=NULL,
                                    "PrecursorCharge"="PrecursorCharge",
                                    "PrecursorMZ"= "FG.PrecMz",
                                    "PrecursorRT"="EG.MeanApexRT",
                                    "PrecursorScore"="EG.Qvalue",
                                    "MS2IntensityAggregated"="FG.Quantity")

.SpectronautFragmentDefsMapping <-c("FragmentIonType"="F.FrgType",
                               "FragmentCharge"="F.FrgZ",
                               "FragmentIntensity"="F.PeakArea")

#'
#' library(readr)
#' res <- read_tsv(file="d:/projects/p2069/data/spectronaut/20160523_130303_20160502_Report.tsv")
#' res %>% glimpse
#' colnames(res)
