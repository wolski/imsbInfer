library(readr)

data2 <- read_tsv("inst/extdata/example.tsv.gz",col_names = TRUE)
dim(data2)

data3 <- prepareOpenSwathData(data)

colnames(data3)

.ProteinDefs <- c("ProteinName",
                  "StrippedSequence","IsotopeLabelType")

.PrecursorDefs <- c("Filename",
                    "StrippedSequence",
                    "ModifiedSequence",
                    "PrecursorCharge",
                    "PrecursorMZ",
                    "PrecursorRT",
                    "PrecursorScore")

.FragmentDefs <-c("Filename",
                  "ModifiedSequence",
                  "PrecursorCharge",
                  "FragmentIonType",
                  "FragmentCharge",
                  "FragmentIntensity")

protein <- data.frame(unique(data[, .ProteinDefs]),stringsAsFactors = FALSE)
precursor <- data.frame(unique(data[, .PrecursorDefs]),stringsAsFactors = FALSE)
transition <- data.frame(unique(data[, .FragmentDefs]),stringsAsFactors = FALSE)
