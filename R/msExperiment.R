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
#' data <- read_tsv("inst/extdata/example.tsv.gz",col_names = TRUE)
#' dim(data)
#' data <- prepareOpenSwathData(data)
#' 
#' huhu <- msTransitionExperiment()
#' huhu$setData(data)
#' intTrans <- huhu$getFragmentIntensities()
#' colnames(intTrans)
#' mypairs(intTrans[,5:ncol(intTrans)],log="xy")
#' precRT <- huhu$getPrecursorRT()
#' library(quantable)

#' mypairs(precRT[,3:ncol(precRT)])
#' precMZ <- huhu$getPrecursorMZ()
#' mypairs(precMZ[,3:ncol(precMZ)])
#' altmanbland(precRT[,3],precRT[,4])
#' precScore <- huhu$getPrecursorScore()
#' mypairs(precScore[,3:ncol(precScore)],log="xy")
#' colnames(huhu$precursor)
#' colnames(huhu$peptide)
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' dim(huhu$peptide)
#' length(unique(huhu$peptide$StrippedSequence))
#' length(unique(huhu$peptide$ProteinName))
#' pep <- huhu$peptide
#' pep<-pep[order(pep$StrippedSequence),]
#' 
#' idx <-which(duplicated(pep$StrippedSequence))
#' pep[c(idx[1]-1,idx[1]),]
#' dim(unique(huhu$peptide[,1:2]))
#' 
#' test <- (huhu$getPrecursorIntensity())
#' huhu$getGlobalFDR()
#' decs<-huhu$getDecoy()

msTransitionExperiment <- 
  setRefClass("msTransitionExperiment",
              fields = list( .data="list", isotopeLabelType = "character",
                             
                             peptide = function(x){
                               if(missing(x)){
                                 return(.data$peptide)
                               }
                               .data$peptide <<- data.frame(unique(x[, .PeptideDefs]),stringsAsFactors = FALSE)
                             },
                             precursor = function(x){
                               if(missing(x)){
                                 return(.data$precursor)
                               }
                               .data$precursor <<-data.frame(unique(x[, .PrecursorDefs]),stringsAsFactors = FALSE)
                             },
                             transition = function(x){
                               if(missing(x)){
                                 return(.data$transition)
                               }
                               .data$transition <<-data.frame(unique(x[, .FragmentDefs]),stringsAsFactors = FALSE)
                             }
              ),
              methods = list(
                initialize = function(...) {
                  require(reshape2)
                  .data$peptide<<-NULL
                  .data$precursor<<-NULL
                  .data$transition<<-NULL
                },
                
                setData = function(data, IsotopeLabelType = "L"){
                  'set the data'
                  data<-data[data$IsotopeLabelType=="L",]
                  isotopeLabelType<<-"L"
                  peptide <<- data
                  precursor <<- data
                  transition <<- data
                },
                
                getFragmentIntensities = function() {
                  'matrix with transitions intensities'
                  transInt <- dcast(transition,
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


