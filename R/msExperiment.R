.PrecursorDefs <- c("Filename",
                    "ProteinName",
                    "Decoy",
                    "StrippedSequence",
                    "ModifiedSequence",
                    "IsotopeLabelType",
                    "PrecursorCharge",
                    "PrecursorMZ",
                    "PrecursorRT",
                    "PrecursorScore"
)

.PeptideDefs <- c("ProteinName",
                  "StrippedSequence",
                  "Decoy")

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

#protein <- data.frame(unique(data[, .ProteinDefs]),stringsAsFactors = FALSE)
#precursor <- data.frame(unique(data[, .PrecursorDefs]),stringsAsFactors = FALSE)
#transition <- data.frame(unique(data[, .FragmentDefs]),stringsAsFactors = FALSE)


sumtop <- function( x , top=3 ){
  if(nrow(x) > top){
    topN = min(nrow(x),top)
    medrow <- apply(x, 1 , median)
    ord<-order(medrow, decreasing = TRUE)[1:topN]
    medrow[ord]
    x<-x[ord,]
  }
  return(apply(x,2,sum,na.rm=TRUE))
}


#' Holds Transition level data
#' 
#' #@field transitiondata transtion data table
#' #@field precdata precursor data table
#' #@field columnsAll required columns for transitiondata
#' #@field columnsPrecursor required columns for precursor
#' @import methods
#' @export msTransitionExperiment
#' @exportClass msTransitionExperiment
#' @examples 
#' rm(list=ls())
#' library(imsbInfer2)
#' library(readr)
#' library(dplyr)
#' data <- read_tsv("inst/extdata/example.tsv.gz",col_names = TRUE)
#' data <- prepareOpenSwathData(data)
#' my_mqdb <- src_sqlite(file.path(".","my_mqdb.sqlite3"), create=TRUE)
#' tmp<-copy_to(my_mqdb, data, "longformat", temporary = FALSE)
#' 
#' huhu <- msTransitionExperiment()
#' huhu$setData(data)
#' intTrans <- huhu$getFragmentIntensities()
#' colnames(intTrans)
#' pairs(intTrans[,5:ncol(intTrans)],log="xy",pch=".")
#' precRT <- huhu$getPrecursorRT()
#' library(quantable)
#' precMZ <- huhu$getPrecursorMZ()
#' precScore <- huhu$getPrecursorScore()
#' mypairs(precScore[,3:ncol(precScore)],log="xy")
#' colnames(huhu$precursor)
#' colnames(huhu$peptide)
#'
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' length(unique(huhu$peptide$StrippedSequence))
#' length(unique(huhu$peptide$ProteinName))
#' 
#' length(idx)
#' pep[c(idx[1]-1,idx[1]),]
#' dim(unique(huhu$peptide[,2:3]))
#' dim(huhu$precursor)
#' dim(merge(huhu$peptide[,c("StrippedSequence","Decoy")], huhu$precursor ))
#' xx <-merge(huhu$peptide[,c("StrippedSequence","Decoy")], huhu$precursor )
#' xx[xx$Decoy ]
#' dim(huhu$precursor)
#' 
#' test <- (huhu$getPrecursorIntensity())
#' huhu$getGlobalFDR()
#' decs<-huhu$getDecoy()
#' 
msTransitionExperiment <- setRefClass("msTransitionExperiment",
              fields = list( .data="src_sqlite",
                             isotopeLabelType = "character", 
                             name="character",
                             path="character",
                             removeDecoy='logical',
                             
                             peptide = function(x){
                               if(missing(x)){
                                 pepcols <- paste(.PeptideDefs, collapse=", ")
                                 query <- c("Select", pepcols, ", count(*) as OrigFreq from LongFormat group by ", pepcols)
                                 query <-paste(tt,collapse=" ")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             },
                             precursor = function(x){
#                                if(missing(x)){
#                                  if( removeDecoy ){
#                                    .precursor <- merge(unique(peptide[,c("StrippedSequence","Decoy")]), .data$precursor )
#                                    .precursor <- .precursor[.precursor$Decoy == 0]
#                                    return(.data$precursor)
#                                  }
#                                  return(.data$precursor)
#                                }
#                                #.data$precursor <<-data.frame(unique(x[, .PrecursorDefs]),stringsAsFactors = FALSE)
                             },
                            transition = function(x){
#                                if(missing(x)){
#                                  if( removeDecoy ){
#                                    .transition <- merge(unique(peptide[,c("StrippedSequence","Decoy")]),.data$transition )
#                                    .transition <- .transition[.transition$Decoy == 0]
#                                    return(.transition)
#                                  }
#                                  return(.data$transition)
#                                }
                               #.data$transition <<-data.frame(unique(x[, .FragmentDefs]),stringsAsFactors = FALSE)
                             }
              ),
              methods = list(
                initialize = function(...) {
                  require(reshape2)
                  isotopeLabelType <<- "L"
                  name <<- uuid::UUIDgenerate()
                  .data <- src_sqlite(file.path(path , name), create=TRUE)
                },
                finilize = function(...){
                  dbDisconnect(.data$con)
                }
                setData = function(data, IsotopeLabelType = "L"){
                  'set the data'
                  library(dplyr)
                  isotopeLabelType<<-"L"
                  data<-data[data$IsotopeLabelType=="L",]
                  tmp<-copy_to(.data, data, "LongFormat", temporary = FALSE)
                  pepcols <- paste(.PeptideDefs, collapse=", ")
                  tt <- c("Select", pepcols, ", count(*) as OrigFreq from LongFormat group by ", pepcols)
                  paste(tt,collapse=" ")  
                },
                noDecoy = function(){
                 removeDecoy <<- TRUE 
                },
                withDecoy = function(){
                 removeDecoy <<- FALSE 
                },
                getFragmentIntensities = function() {
                  'matrix with transitions intensities'
                  
                  transInt <- dcast(.transition,
                                    ModifiedSequence + PrecursorCharge + FragmentIonType + FragmentCharge ~ Filename ,
                                    value.var="FragmentIntensity")
                  return(transInt)
                },
                
                
                getPrecursorRT = function() {
                  'matrix with precursor retention times'
                  transRT = dcast(precursor, ModifiedSequence + PrecursorCharge  ~ Filename , value.var="PrecursorRT")
                },
                
                getPrecursorScore = function() {
                  'matrix with precursor scores'
                  transScore = dcast(precursor, ModifiedSequence + PrecursorCharge ~ Filename , value.var="PrecursorScore")
                  return(transScore)
                },
                getPrecursorMZ = function() {
                  'matrix with precursor mz'
                  transMz = dcast(precursor, ModifiedSequence + PrecursorCharge  ~ Filename , value.var="PrecursorMZ")
                  return(transMz)
                },
                
              
                getPrecursorIntensity=function(FUN=sum){
                  transIntensities <- getFragmentIntensities()
                  key <- paste(transIntensities$ModifiedSequence, ModifiedSequence$z, sep="_")
                  res <- by(transIntensities,INDICES=key)
                  return(res)
                },
                
                getGlobalFDR = function(){
                  'compute FDR for dataset'
                  
                  tmp <- precursor[, c("ProteinName", "Decoy", "PrecursorScore")]
                  tmp <- tmp[tmp$Score < 2,] # not sure how to treat requant values in this context
                  fdr <- (sum(tmp$Decoy) / length(tmp$ProteinName))
                  return(fdr)
                },
                plotSampleFDR = function(log=""){
                  'plot FDR versus Score'
                  tmp <-data.frame(Protein = precdata$Protein,Decoy = precdata$Decoy,Score = precdata$Score)
                  tmp <- tmp[tmp$Score < 2,] # not sure how to treat requant values in this context
                  tmp <- tmp[order(tmp$Score),]
                  plot(tmp$Score,cumsum(tmp$Decoy) / nrow(tmp) * 100 ,type="l",xlab="Score", ylab="FDR",log=log)
                }
                
              ))


