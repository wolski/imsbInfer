.PrecursorDefs <- c("FileName",
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

.PrecursorDefs <- c("FileName",
                    "StrippedSequence",
                    "ModifiedSequence",
                    "PrecursorCharge",
                    "PrecursorMZ",
                    "PrecursorRT",
                    "PrecursorScore")

.FragmentDefs <-c("FileName",
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
#' "src_sqlite" class
#'
#' @name src_sqlite
#' @aliases src_sqlite
#' @family src_sqlite
#'
#' @exportClass src_sqlite
setOldClass("src_sqlite")
#' Holds Transition level data
#' @field name name of experiment
#' @field path where to store data
#' @field peptide transtion data table
#' @field precursor precursor data table
#' @field transition required columns for transitiondata
#' @import methods
#' @import dplyr
#' @import reshape2
#' @import RSQLite
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
#' huhu <- msTransitionExperiment(path=".", name="mydb1.sql")
#' huhu$setData(data)
#' xx <-huhu$name
#' huhu$getFileName()
#' 
#' huhu$finalize()
#' 
#' huhu <- msTransitionExperiment(path=".", name="mydb.sql")
#' intTrans <- huhu$getFragmentIntensities()
#' dim(intTrans)
#' precInt <- huhu$getPrecursorIntensitySum()
#' dim(precInt)
#' colnames(precInt)
#' pairs(precInt[,3:ncol(precInt)], pch='.',log="xy")
#' huhu$noDecoy()
#' colnames(intTrans)
#' intTrans <- huhu$getFragmentIntensities()
#' head(intTrans)
#' xx <- by(intTrans[,4:ncol(intTrans)], intTrans$ModifiedSequence , FUN = function(x){x})
#' length(xx)
#' dim(xx[[1]])
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
#' 
#' test <- (huhu$getPrecursorIntensity())
#' huhu$getGlobalFDR()
#' decs<-huhu$getDecoy()
#'

msTransitionExperiment <-
  setRefClass("msTransitionExperiment",
              fields = list( .data="src_sqlite",
                             isotopeLabelType = "character", 
                             name="character",
                             path="character",
                             .removeDecoy='logical',
                             
                             peptide = function(x){
                               print("in peptide")
                               if(missing(x)){
                                 where <- ""
                                 if(.removeDecoy){
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
                                 if(.removeDecoy){
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
                                 if(.removeDecoy){
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
                                      .removeDecoy=FALSE,
                                      isotopeLabelType="L",
                                      ...) {
                  print(match.call())
                  require(reshape2)
                  
                  name <<- name
                  path <<- path
                  isotopeLabelType <<- isotopeLabelType
                  .removeDecoy <<- .removeDecoy
                  
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
                  (dplyr::collect(dplyr::tbl(.data,dplyr::sql("SELECT * FROM LongFormat"))))
                },
                getFileName = function(){
                  dplyr::collect(dplyr::tbl(.data,dplyr::sql("SELECT DISTINCT FileName FROM LongFormat")))
                },
                noDecoy = function(){
                  .removeDecoy <<- TRUE 
                },
                withDecoy = function(){
                  .removeDecoy <<- FALSE 
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

                getPrecursorIntensitySum=function(){
                  where <- ""
                  if(.removeDecoy){
                    where <- " where Decoy = 0 "
                  }
                  
                  fragmentCols <- paste(.FragmentDefs, collapse=", ")
                  query <- c("Select Filename, ModifiedSequence, PrecursorCharge, Decoy, count(*) as Freq, sum(FragmentIntensity) as Intensity from LongFormat " ,
                             where , " group by  Filename, ModifiedSequence, PrecursorCharge")
                  query <-paste(query,collapse=" ")
                  print(query)
                  tmp <- dbGetQuery(.data$con,query)
                  transPrecInt = dcast(tmp, ModifiedSequence + PrecursorCharge  ~ Filename , value.var="Intensity")
                  return(transPrecInt)
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


