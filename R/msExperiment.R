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
#' 
#' huhu <- msTransitionExperiment()
#' huhu$setData(data)
#' data <- huhu$getData() 
#' x<-huhu$peptide
#' head(x)
#' x <- huhu$
#' dumm <- msTransitionExperiment(path=".", name="mydb.sql")
#' dumm$name
#' dumm$setData(data)
#' huhu <- dumm
#' intTrans <- huhu$getFragmentIntensities()
#' dim(intTrans)
#' huhu$noDecoy()
#' colnames(intTrans)
#' intTrans <- huhu$getFragmentIntensities()
#' dim(intTrans)
#' pairs(intTrans[,5:ncol(intTrans)],log="xy",pch=".")
#' precRT <- huhu$getPrecursorRT()
#' precMZ <- huhu$getPrecursorMZ()
#' precScore <- huhu$getPrecursorScore()
#' pairs(precScore[,3:ncol(precScore)],log="xy",pch=".")
#' colnames(huhu$precursor)
#' colnames(huhu$peptide)
#'
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' huhu$withDecoy()
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' huhu$noDecoy()
#' length(unique(huhu$peptide$StrippedSequence))
#' length(unique(huhu$peptide$ProteinName))
#' xx <-merge(huhu$peptide[,c("StrippedSequence","Decoy")], huhu$precursor )
#' head(xx)
#' 
#' test <- (huhu$getPrecursorIntensity())
#' huhu$getGlobalFDR()
#' decs<-huhu$getDecoy()
#' 
setOldClass("src_sqlite")
msTransitionExperiment <-
  setRefClass("msTransitionExperiment",
              fields = list( .data="src_sqlite",
                             isotopeLabelType = "character", 
                             name="character",
                             path="character",
                             removeDecoy='logical',
                             
                             peptide = function(x){
                               print("in peptide")
                               if(missing(x)){
                                 where <- ""
                                 if(removeDecoy){
                                   where <- " where Decoy = 0 "
                                 }
                                 peptideCols <- paste(.PeptideDefs, collapse=", ")
                                 query <- c("Select", peptideCols, ", count(*) as Freq from LongFormat ", where , " group by ", peptideCols)
                                 query <-paste(query,collapse=" ")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             },
                             precursor = function(x){
                               if(missing(x)){
                                 where <- ""
                                 if(removeDecoy){
                                   where <- " where Decoy = 0 "
                                 }
                                 
                                 precursorCols <- paste(.PrecursorDefs, collapse=", ")
                                 query <- c("Select", precursorCols, ", count(*) as Freq from LongFormat " , where , " group by ", precursorCols)
                                 query <-paste(query,collapse=" ")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             },
                             transition = function(x){
                               if(missing(x)){
                                 where <- ""
                                 if(removeDecoy){
                                   where <- " where Decoy = 0 "
                                 }
                                 
                                 fragmentCols <- paste(.FragmentDefs, collapse=", ")
                                 query <- c("Select", fragmentCols, ", count(*) as Freq from LongFormat " , where , " group by ", fragmentCols)
                                 query <-paste(query,collapse=" ")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             }
              ),
              methods = list(
                initialize = function(name=paste( uuid::UUIDgenerate(), ".sqlite",sep=""),
                                      path=".",
                                      removeDecoy=FALSE,
                                      isotopeLabelType="L",
                                      ...) {
                  print(match.call())
                  require(reshape2)
                  
                  name <<- name
                  path <<- path
                  isotopeLabelType <<- isotopeLabelType
                  removeDecoy <<- removeDecoy
                  
                  dbfile <- file.path(path , name)
                  print(dbfile)
                  .data <<- src_sqlite(dbfile, create=TRUE)
                  print('done')
                },
                finalize  = function(){
                  print("infinalize")
                  require("RSQLite")
                  dbDisconnect(.self$.data$con)
                  print("disconnected")
                },
                setData = function(data, IsotopeLabelType = "L"){
                  'set the data'
                  library(dplyr)
                  isotopeLabelType<<-"L"
                  data<-data[data$IsotopeLabelType=="L",]
                  tmp <- copy_to(.data, data, "LongFormat", temporary = FALSE)
                  dbListTables(.data$con)
                },
                getData = function(){
                  print(dbListTables(.data$con))
                  
                  as.data.frame(dplyr::collect(dplyr::tbl(.data,dplyr::sql("SELECT * FROM LongFormat"))))
                },
                noDecoy = function(){
                  removeDecoy <<- TRUE 
                },
                withDecoy = function(){
                  removeDecoy <<- FALSE 
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
                  "!!!don't use. work in progress!!!"
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


