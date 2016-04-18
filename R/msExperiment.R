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
#' @export msTransitions
#' @exportClass msTransitions
#' @examples 
#' rm(list=ls())
#' library(imsbInfer2)
#' library(readr)
#' 
#' data <- read_tsv(file,col_names = TRUE)
#' data <- prepareOpenSwathData(data)
#' huhu <- msTransitions()
#' huhu$setData(prepOS)
#' huhu$getFDR()
#' huhu$plotFDR()
#' decs<-huhu$getDecoy()

msTransitionExperiment <- setRefClass("msTransitionsExperiment",
                                      fields = list( 
                                        protein = function(x){
                                          data.frame(unique(x[, .ProteinDefs]),stringsAsFactors = FALSE)
                                        },
                                        precursor = function(x){
                                          data.frame(unique(x[, .PrecursorDefs]),stringsAsFactors = FALSE)
                                        },
                                        transition = function(x){
                                          data.frame(unique(data[, .FragmentDefs]),stringsAsFactors = FALSE)
                                        }
                                      ),
                                      methods = list(
                                        initialize = function(...) {
                                          require(reshape2)
                                        },
                                        setData = function(data){
                                          'set the data'
                                          protein <<- data
                                          precursor <<- data
                                          transition <<- data
                                        },
                                        
                                        getFragmentIntensities = function() {
                                          'matrix with transitions intensities'
                                          transIntensities <- dcast(transition,
                                                                    ModifiedSequence + PrecursorCharge + FragmentIonType + FragmentCharge ~ Filename ,
                                                                    value.var="FragmentIntensity")
                                          return(transIntensities)
                                        },
                                        
                                        getPrecursorIntensity=function(FUN=sum){
                                          transIntensities <- getFragmentIntensities()
                                          key <- paste(transIntensities$ModifiedSequence, ModifiedSequence$z, sep="_")
                                          res <- by(transIntensities,INDICES=key,FUN=FUN)
                                          do.call("rbind",res)
                                        },
                                        
                                        getPrecursorRT = function() {
                                          'matrix with precursor retention times'
                                          transRT = dcast(precdata, ModifiedSequence + PrecursorCharge  ~ Filename , value.var="PrecursorRT")
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
                                        getDecoy= function(){
                                          'matrix with decoy information'
                                          transDecoy = dcast(precursor, ModifiedSequence + PrecursorCharge  ~ Filename , value.var="Decoy")
                                          return(transDecoy)
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


