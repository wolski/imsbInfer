---
  title: "Play along"
author: "Witold Wolski"
date: "11 March 2016"
output: html_document
---


```{r}

split2long <- function(far){
  far <- data2
  #far = feature_alignment_requant
  apa = as.character(far$aggr_Peak_Area)
  afa = as.character(far$aggr_Fragment_Annotation)
  # split transition intensities
  transints = lapply(apa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # split transition names
  transids = lapply(afa,function(x){unlist(strsplit(x,";",fixed=TRUE))})
  # prepare output
  lx = length(transids)
  stopifnot(lx == nrow(far))

  far<- far[, !(names(far) %in% c("aggr_Peak_Area","aggr_Fragment_Annotation"))]
  message("prepared dataframe")
  # extend
  lengths <-sapply(transids,length)
  idx <- rep(1:lx, lengths)
  range(idx)
  far <- far[idx, ]
  message("and finally nrow " , nrow(far), " ncol ", ncol(far))
  far <- data.frame(far, aggr_Peak_Area = as.numeric(unlist(transints)),
                    aggr_Fragment_Annotation =  as.character(unlist( transids) ) )
  return(far)
}

prepareDF <- function (df)
{
  colnames(df) <- gsub("m/z","mz",colnames(df))
  required = c("transition_group_id", "align_origfilename","decoy",
               "RT", "mz", "Intensity", "ProteinName", "m_score", "aggr_Fragment_Annotation",
               "aggr_Peak_Area" )
  x = match(tolower(required), tolower(colnames(df)))
  stopifnot(required == colnames(df)[x])
  df = df[, x]
  df <- df[order(df$transition_group_id), ]
  return(df)
}

```

```{r}
rm(list=ls())
gc()
library(imsbInfer)
library(readr)

SpecLib = "data/E1509221443_feature_alignment_requant.tsv"

data <- read_tsv(SpecLib,col_names = TRUE)
#tmp <-read.csv(SpecLib, stringsAsFactors = F,sep="\t",header=T)


data2 <- prepareDF(data)
data2$align_origfilename <- gsub("_with_dscore_filtered.csv","", basename(data2$align_origfilename))

plot(table(data2$decoy))


data3 <- split2long(data2)



transIntensities = dcast(data3, transition_group_id + aggr_Fragment_Annotation ~ align_origfilename , value.var="aggr_Peak_Area")


transIntensities$aggr_Fragment_Annotation <- gsub("^DECOY_","DECOY",transIntensities$aggr_Fragment_Annotation )
tmp <-split2table(as.character(transIntensities$aggr_Fragment_Annotation) , split = "_")
colnames(tmp)<-c("frag_id", "ion_type","fragCharge","pepSequence","pepCharge")
head(transIntensities)


mscore = dcast(data2, transition_group_id  ~ align_origfilename , value.var="m_score")
rt = dcast(data2, transition_group_id ~ align_origfilename , value.var="RT")
mz = dcast(data2, transition_group_id ~ align_origfilename , value.var="mz")
pepIntensities = dcast(data2,  transition_group_id ~ align_origfilename , value.var="Intensity")

head(data3)
split2table(data)

```
