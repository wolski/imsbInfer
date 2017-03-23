
.PeptideDefs <- c("ProteinName",
                  "StrippedSequence",
                  "Decoy")

.PrecursorDefs <- c("FileName",
                    "StrippedSequence",
                    "ModifiedSequence",
                    "PrecursorCharge",
                    "PrecursorMZ",
                    "PrecursorRT",
                    "LabelType",
                    "Decoy",
                    "PrecursorScore")

.PrecursorIntensity <- c("FileName",
                         "StrippedSequence",
                         "ModifiedSequence",
                         "PrecursorCharge",
                         "MS1Intensity", # LFQ
                         "MS2IntensityAggregated"
)

#Dia stuff
.FragmentDefs <-c("FileName",
                  "ModifiedSequence",
                  "PrecursorCharge",
                  "FragmentIonType",
                  "FragmentCharge",
                  "FragmentInterference",
                  "FragmentIntensity")


.RequiredColumns <- base::unique(c(.PeptideDefs, .PrecursorDefs,.PrecursorIntensity, .FragmentDefs ))
.RequiredColumns

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
#' @field labelType Stores only one label type
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
#' data <- read_tsv(file.path(path.package("imsbInfer2"),"extdata/example.tsv.gz"),col_names = TRUE)
#' data <- prepareOpenSwathData(data)
#' dim(data)
#' 
#' huhu <- msTransitionExperiment(path=".", name="mydb1.sql")
#' head(data)
#' huhu$setData(data)
#' 
#' huhu$name
#' huhu$getFileName()
#' 
#' 
#' intTrans <- huhu$getFragmentIntensities()
#' dim(intTrans)
#' head(intTrans)
#' precInt <- huhu$getPrecursorIntensitySum()
#' dim(precInt)
#' colnames(precInt)
#' pairs(precInt[,3:ncol(precInt)], pch='.',log="xy")
#' 
#' precIntExt <- huhu$getPrecursorIntensity()
#' dim(precIntExt)
#' head(precInt) == head(precIntExt)
#' 
#' huhu$noDecoy()
#' intTrans <- huhu$getFragmentIntensities()
#' xx <- by(intTrans[,5:ncol(intTrans)], intTrans$ModifiedSequence , FUN = function(x){x})
#' length(xx)
#' xx[[1]]
#' pairs(intTrans[,5:ncol(intTrans)],log="xy",pch=".")
#' precRT <- huhu$getPrecursorRT()
#' precMZ <- huhu$getPrecursorMZ()
#' precScore <- huhu$getPrecursorScore()
#' pairs(precScore[,3:ncol(precScore)],log="xy",pch=".")
#' colnames(huhu$precursor)
#' head(huhu$peptide$Freq)
#'
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' huhu$withDecoy()
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' huhu$noDecoy()
#' length(unique(huhu$peptide$StrippedSequence))
#' length(unique(huhu$peptide$ProteinName))
#' xx <-merge(huhu$peptide[,c("StrippedSequence","Decoy")], huhu$precursor )
#' head(xx)
#' dim(xx)
#' huhu$withDecoy()
#' huhu$getGlobalFDR()
#' huhu$plotFalsePostives()
#'
msTransitionExperiment <-
  setRefClass("msTransitionExperiment",
              fields = list( .data="src_sqlite",
                             labelType = "character", 
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
                                      labelType="L",
                                      ...) {
                  print(match.call())
                  require(reshape2)
                  
                  name <<- name
                  path <<- path
                  labelType <<- labelType
                  .removeDecoy <<- .removeDecoy
                  
                  dbfile <- file.path(path , name)
                  
                  .data <<- src_sqlite(dbfile, create=TRUE)
                  print('done')
                },
                
                finalize  = function(){
                  print("infinalize")
                  require("RSQLite")
                  dbDisconnect(.self$.data$con)
                  print("disconnected")
                },
                
                setData = function(data, labelType = "L"){
                  'set the data (drops all data if already existing)'
                  stopifnot(.RequiredColumns %in% colnames(data))
                  library(dplyr)
                  labelType <<- labelType
                  data<-data[data$LabelType==labelType,]
                  print(dim(data))
                  if(db_has_table(.data$con, "LongFormat")){
                    dplyr::db_drop_table(.data$con,  "LongFormat", force = FALSE)
                  }
                  tmp <- dplyr::copy_to(.data, data, "LongFormat", temporary = FALSE)
                  print("data inserted")
                  dbListTables(.data$con)
                },
                getData = function(){
                  'get all the data in long format'
                  print(dbListTables(.data$con))
                  (dplyr::collect(dplyr::tbl(.data,dplyr::sql("SELECT * FROM LongFormat"))))
                },
                getFileName = function(){
                  'get filenames'
                  dplyr::collect(dplyr::tbl(.data,dplyr::sql("SELECT DISTINCT FileName FROM LongFormat")))
                },
                
                noDecoy = function(){
                  'work in decoy free mode'
                  .removeDecoy <<- TRUE 
                },
                withDecoy = function(){
                  'work with decoys'
                  .removeDecoy <<- FALSE 
                },
                getFragmentIntensities = function() {
                  'matrix with transitions intensities'
                  
                  transInt <- dcast(transition,
                                    ModifiedSequence + PrecursorCharge + FragmentIonType + FragmentCharge ~ FileName ,
                                    value.var="FragmentIntensity")
                  return(transInt)
                },
                
                
                getPrecursorRT = function() {
                  'matrix with precursor retention times'
                  transRT = dcast(precursor, ModifiedSequence + PrecursorCharge  ~ FileName , value.var="PrecursorRT")
                },
                
                getPrecursorScore = function() {
                  'matrix with precursor scores'
                  transScore = dcast(precursor, ModifiedSequence + PrecursorCharge ~ FileName , value.var="PrecursorScore")
                  return(transScore)
                },
                getPrecursorMZ = function() {
                  'matrix with precursor mz'
                  transMz = dcast(precursor, ModifiedSequence + PrecursorCharge  ~ FileName , value.var="PrecursorMZ")
                  return(transMz)
                },
                
                getPrecursorIntensitySum=function(){
                  'Selects Fragment intensities and aggregates them (sum)'
                  
                  where <- NULL 
                  if(.removeDecoy){
                    where <- " where Decoy = 0 "
                  }
                  fragmentCols <- paste(.FragmentDefs, collapse=", ")
                  query <- c("Select FileName, ModifiedSequence, PrecursorCharge, Decoy, count(*) as Freq, sum(FragmentIntensity) as Intensity from LongFormat " ,
                             where , " group by  FileName, ModifiedSequence, PrecursorCharge")
                  query <-paste(query,collapse=" ")
                  print(query)
                  tmp <- dbGetQuery(.data$con,query)
                  transPrecInt = dcast(tmp, ModifiedSequence + PrecursorCharge  ~ FileName , value.var="Intensity")
                  return(transPrecInt)
                },
                getPrecursorIntensity=function(){
                  'Selects Fragment intensities given by external software'
                  
                  where <- NULL 
                  if(.removeDecoy){
                    where <- " where Decoy = 0 "
                  }
                  fragmentCols <- paste(.FragmentDefs, collapse=", ")
                  query <- c("Select Distinct FileName, ModifiedSequence, PrecursorCharge, Decoy, MS2IntensityAggregated as Intensity from LongFormat " ,
                             where )
                  query <-paste(query,collapse=" ")
                  print(query)
                  tmp <- dbGetQuery(.data$con,query)
                  transPrecInt = dcast(tmp, ModifiedSequence + PrecursorCharge  ~ FileName , value.var="Intensity")
                  return(transPrecInt)
                },
                getGlobalFDR = function(){
                  'compute FDR for dataset'
                  tmp <- precursor[, c("Decoy", "PrecursorScore")]
                  tmp <- tmp[tmp$PrecursorScore < 2,] # not sure how to treat requant values in this context
                  fdr <- (sum(tmp$Decoy) / nrow(tmp))
                  return(fdr)
                },
                
                plotFalsePostives = function(log=""){
                  'plot FP versus Score'
                  tmp <- precursor[, c("Decoy", "PrecursorScore")]
                  tmp <- tmp[tmp$PrecursorScore < 2,] # not sure how to treat requant values in this context
                  tmp <- tmp[order(tmp$PrecursorScore),]
                  plot(tmp$PrecursorScore,cumsum(tmp$Decoy) / nrow(tmp) * 100 ,type="l",xlab="Score", ylab="FP",log=log)
                }
              )
  )


