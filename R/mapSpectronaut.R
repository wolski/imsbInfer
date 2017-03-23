.SpectronautPrecursorDefsMapping <- list("FileName"="Run",
                                    "ProteinName"="ProteinName",
                                    "Decoy"="EG.IsDecoy",
                                    "StrippedSequence"="EG.StrippedSequence",
                                    "ModifiedSequence"="PeptideSequence",
                                    "PrecursorCharge"="PrecursorCharge",
                                    "PrecursorMZ"= "FG.PrecMz",
                                    "PrecursorRT"="EG.MeanApexRT",
                                    "PrecursorScore"="EG.Qvalue",
                                    "MS2IntensityAggregated"="FG.Quantity")

.SpectronautFragmentDefsMapping <- c("FragmentIonType"="F.FrgType",
                               "FragmentCharge"="F.FrgZ",
                               "FragmentIntensity"="F.PeakArea")

names(.OpenMSPrecursorDefsMapping) == names(.SpectronautPrecursorDefsMapping)
names(.SpectronautFragmentDefsMapping) == names(.OpenMSFragmentDefsMapping)

#'
#' library(readr)
#' res <- read_tsv(file="d:/projects/p2069/data/spectronaut/20160523_130303_20160502_Report.tsv")
#' res %>% glimpse
#' colnames(res)
#' res <- read_tsv(file="d:/projects/p1342_MarianaPardo_DIA/20170309_171629_p1342_marianaPardo_Report_TRANSITION.xls")
#'
#' res <- read_tsv(file="d:/projects/p2069_WolfgangFaigle_Citr_BR/dataQuantitative/spectronaut/20160523_130303_20160502_Report.tsv")
#' .SpectronautPrecursorDefsMapping %in% colnames(res)
#' .SpectronautFragmentDefsMapping %in% colnames(res)
#' colnames(res)[!colnames(res) %in% c(.SpectronautPrecursorDefsMapping,.SpectronautFragmentDefsMapping )]
