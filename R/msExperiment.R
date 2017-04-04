
.PeptideDefs <- c("ProteinName",
                  "StrippedSequence",
                  "Decoy")

.PrecursorID <- c(  "TransitionGroupID",
                    "StrippedSequence",
                    "ModifiedSequence",
                    "PrecursorCharge",
                    "Decoy")

.PrecursorDefs <- c("FileName",
                    "TransitionGroupID",
                    "StrippedSequence",
                    "ModifiedSequence",
                    "PrecursorCharge",
                    "PrecursorMZ",
                    "PrecursorRT",
                    "LabelType",
                    "Decoy",
                    "PrecursorScore")

.PrecursorIntensity <- c("FileName",
                         "TransitionGroupID",
                         "StrippedSequence",
                         "ModifiedSequence",
                         "PrecursorCharge",
                         "MS1Intensity", # LFQ
                         "MS2IntensityAggregated"
)

#Dia stuff
.FragmentDefs <-c("FileName",
                  "TransitionGroupID",
                  "FragmentID",
                  "ModifiedSequence",
                  "PrecursorCharge",
                  "FragmentIonType",
                  "FragmentCharge",
                  "FragmentInterference",
                  "FragmentIntensity")


.RequiredColumns <- base::unique(c(.PeptideDefs, .PrecursorDefs,.PrecursorIntensity, .FragmentDefs ))
.PrecursorColumns <- base::unique(c(.PeptideDefs, .PrecursorDefs,.PrecursorIntensity))

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
#' 
#' huhu <- msTransitionExperiment(path=".", name="mydb2.sql")
#' huhu$setData(data)
#' tmp <- huhu$getData()
#' stopifnot(dim(tmp) == dim(data))
#' huhu$getFileName()
#' 
#' intTrans <- huhu$getFragmentIntensities()
#' dim(intTrans)
#' head(intTrans)
#' 
#' precIntExt <- huhu$getPrecursorIntensity()
#' head(precIntExt)
#' dim(precIntExt)
#' huhu$noDecoy()
#' 
#' precRT <- huhu$getPrecursorRT()
#' precMZ <- huhu$getPrecursorMZ()
#' precScore <- huhu$getPrecursorScore()
#' head(precScore)
#' pairs(precScore[,6:ncol(precScore)],log="xy",pch=".")
#' dim(huhu$precursor)
#' dim(huhu$precursorId)
#' table(huhu$peptide$Freq)
#'
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' huhu$withDecoy()
#' head(huhu$peptide[huhu$peptide$Decoy==1,])
#' dim(huhu$precursor)
#' dim(huhu$precursorId)
#' huhu$noDecoy()
#' length(unique(huhu$peptide$StrippedSequence))
#' length(unique(huhu$peptide$ProteinName))
#' xx <-merge(huhu$peptide[,c("StrippedSequence","Decoy")], huhu$precursor )
#' head(xx)
#' huhu$withDecoy()
#' huhu$getGlobalFDR()
#' huhu$plotFalsePostives()
msTransitionExperiment <-
  setRefClass("msTransitionExperiment",
              fields = list( .data="src_sqlite",
                             labelType = "character", 
                             name="character",
                             path="character",
                             .removeDecoy='logical',
                             
                             peptide = function(x){
                               message("in peptide")
                               if(missing(x)){
                                 where <- ""
                                 if(.removeDecoy){
                                   where <- " where Decoy = 0 "
                                 }
                                 peptideCols <- paste(.PeptideDefs, collapse=", ")
                                 query <- c("Select", peptideCols, ", count(*) as Freq from PrecursorInfo ", where , " group by ", peptideCols)
                                 query <-paste(query,collapse=" ")
                                 message(query)
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             },
                             precursorId = function(x){
                               message("in precursor")
                               if(missing(x)){
                                 where <- ""
                                 if(.removeDecoy){
                                   where <- " where Decoy = 0 "
                                 }
                                 precursorIDCols <- paste(.PrecursorID, collapse=", ")
                                 query <- c("Select", precursorIDCols, ", count(*) as Freq from PrecursorInfo ", where , " group by ", precursorIDCols)
                                 query <-paste(query,collapse=" ")
                                 message(query)
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             },
                             precursor = function(x){
                               if(missing(x)){
                                 where <- ""
                                 if(.removeDecoy){
                                   where <- " where Decoy = 0 "
                                 }
                                 
                                 precursorCols <- paste(.PrecursorColumns, collapse=", ")
                                 query <- c("Select", precursorCols, ", count(*) as Freq from PrecursorInfo " , where , " group by ", precursorCols)
                                 query <-paste(query,collapse=" ")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             },
                             transition = function(x){
                               if(missing(x)){
                                 fragmentCols <- paste(.FragmentDefs, collapse=", ")
                                 query <- c("Select", fragmentCols, ", count(*) as Freq from FragmentInfo group by ", fragmentCols)
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
                  message(match.call())
                  require(reshape2)
                  
                  name <<- name
                  path <<- path
                  labelType <<- labelType
                  .removeDecoy <<- .removeDecoy
                  
                  dbfile <- file.path(path , name)
                  
                  .data <<- src_sqlite(dbfile, create=TRUE)
                  message('done')
                },
                finalize  = function(){
                  message("infinalize")
                  require("RSQLite")
                  dbDisconnect(.self$.data$con)
                  message("disconnected")
                },
                setData = function(data, labelType = "L"){
                  'set the data (drops all data if already existing)'
                  stopifnot(.RequiredColumns %in% colnames(data))
                  library(dplyr)
                  labelType <<- labelType
                  data<-data[data$LabelType==labelType,]
                  message(nrow(data), " ", ncol(data))
                  
                  precursorInfo <- unique(data[, .PrecursorColumns])
                  fragmentInfo <- unique(data[, .FragmentDefs])
                  
                  if(db_has_table(.data$con, "FragmentInfo")){
                    dplyr::db_drop_table(.data$con,  "FragmentInfo", force = FALSE)
                  }
                  tmp <- dplyr::copy_to(.data, fragmentInfo, "FragmentInfo", temporary = FALSE)

                  if(dplyr::db_has_table(.data$con, "PrecursorInfo")){
                    dplyr::db_drop_table(.data$con,  "PrecursorInfo", force = FALSE)
                  }
                  tmp <- dplyr::copy_to(.data, precursorInfo, "PrecursorInfo", temporary = FALSE)
                  
                  res <- dbListTables(.data$con)
                  message(res)
                  (res)
                },
                getData = function(){
                  'get all the data in long format'
                  PrecInfo <- DBI::dbGetQuery(.data$con,"Select * from PrecursorInfo")
                  FragmentInfo <- dbGetQuery(.data$con,"Select * from FragmentInfo")
                  merge(PrecInfo,FragmentInfo )
                },
                getFileName = function(){
                  'get filenames'
                  dplyr::collect(dplyr::tbl(.data,dplyr::sql("SELECT DISTINCT FileName FROM PrecursorInfo")))
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
                                    TransitionGroupID + FragmentID + ModifiedSequence + PrecursorCharge + FragmentIonType + FragmentCharge ~ FileName ,
                                    value.var="FragmentIntensity")
                  return(transInt)
                },
                dcastPrecursor = function(x, value.var="Intensity"){
                  message("casting using value.var : ", value.var)
                  dcast(x, TransitionGroupID + ModifiedSequence + PrecursorCharge + StrippedSequence + Decoy ~ FileName , value.var=value.var, drop=TRUE)
                },
                acastPrecursor = function(x, value.var="Intensity"){
                  message("casting using value.var : ", value.var)
                  acast(x, TransitionGroupID + ModifiedSequence + PrecursorCharge + StrippedSequence + Decoy ~ FileName , value.var=value.var, drop=TRUE)
                },
                
                getPrecursorRT = function() {
                  'matrix with precursor retention times'
                  transRT = dcastPrecursor(precursor, value.var="PrecursorRT")
                },
                
                getPrecursorScore = function() {
                  'matrix with precursor scores'
                  transScore = dcastPrecursor(precursor,  value.var="PrecursorScore")
                  return(transScore)
                },
                getPrecursorMZ = function() {
                  'matrix with precursor mz'
                  transMz = dcastPrecursor(precursor,  value.var="PrecursorMZ")
                  return(transMz)
                },
                # getPrecursorIntensitySum=function(){
                #   'Selects Fragment intensities and aggregates them (sum)'
                #   where <- NULL 
                #   if(.removeDecoy){
                #     where <- " where Decoy = 0 "
                #   }
                #   fragmentCols <- paste(.FragmentDefs, collapse=", ")
                #   query <- c("Select TransitionGroupID, FileName, ModifiedSequence, StrippedSequence, PrecursorCharge, Decoy, count(*) as Freq, sum(FragmentIntensity) as Intensity from FragmentInfo " ,
                #              where , 
                #              " group by TransitionGroupID, FileName, ModifiedSequence, PrecursorCharge")
                #   query <-paste(query,collapse=" ")
                #   message(query)
                #   tmp <- dbGetQuery(.data$con,query)
                #   #tmp <- dcastPrecursor(tmp, value.var = "Intensity")
                #   return(tmp)
                # },
                getPrecursorIntensity = function(){
                  #TODO(WEW): check why different result than getPrecursorIntensitySum
                  'DO NOT USE!!! Selects Fragment intensities given by external software, '
                  where <- NULL
                  if(.removeDecoy){
                    where <- " where Decoy = 0 "
                  }
                  fragmentCols <- paste(.PrecursorIntensity , collapse=", ")
                  query <- c("Select",fragmentCols  ,  " from PrecursorInfo " , where, "group by", fragmentCols  )
                  query <-paste(query,collapse=" ")
                  message(query)
                  tmp <- dbGetQuery(.data$con,query)
                  plyr::rename(tmp, c("MS2IntensityAggregated" = "Intensity"))
                  return(tmp)
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


