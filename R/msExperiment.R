.formalAllDef <-c("z","Protein", "RT",  "ModifiedSequence",
                 "Sequence", "Filename", "Decoy",
                 "mz", "Score", "IonType",
                 "FragZ", "Intensity")

.formalPrecDef <- c("z", "RT","Protein",  "ModifiedSequence",
                   "Sequence", "Filename", "Decoy",
                   "mz", "Score")


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
#' @field transitiondata transtion data table
#' @field precdata precursor data table
#' @field columnsAll required columns for transitiondata
#' @field columnsPrecursor required columns for precursor
#' @import methods
#' @export msTransitions
#' @exportClass msTransitions
#' @examples 
#' rm(list=ls())
#' library(imsbInfer2)
#' library(readr)
#' huhu$columnsAll
#' huhu$columnsPrecursor
#' huhu$transitiondata
#' file = "C:/Users/wewol/Google Drive/tissuecomparison/OpenSWATH/BAT_19strains/data/E1603291025_feature_alignment_requant.tsv.gz"
#' data <- read_tsv(file,col_names = TRUE)
#' data <- prepareOpenSwathData(data)
#' head(data)
#' 
#' prepOS <- data
#' save(prepOS, file="inst/temp/prepOS.Rd")
#' load(file="inst/temp/prepOS.Rd")
#' huhu <- msTransitions()
#' huhu$setData(prepOS)
#' huhu$getFDR()
#' huhu$plotFDR()
#' decs<-huhu$getDecoy()
#' head(decs)
#' head(decs)
#' xx <- apply(decs[,4:ncol(decs)] , 1, sum)
#' table(xx)[2]/table(xx)[1] * 100
#' 
#' intensity <- huhu$getTransIntensity()
#' dim(intensity)
#' rt <-huhu$getRT()

"c:/Users/Google "
msTransitions <- setRefClass("msTransitions",
                             fields = list( transitiondata = "data.frame",
                                            precdata = "data.frame",
                                            columnsAll="character",
                                            columnsPrecursor="character"),
                             methods = list(
                               initialize = function(...) {
                                 require(reshape2)
                                 columnsAll <<- .formalAllDef
                                 columnsPrecursor <<- .formalPrecDef
                               },
                               
                               .validateColumns = function(colnames){
                                 if(sum(columnsAll %in% colnames) != length(columnsAll)){
                                   message("missing columns : ")
                                   message(columnsAll[!(columnsAll %in% colnames)])
                                   return(FALSE)
                                   
                                 }
                                 return(TRUE)
                               },
                               
                               setData = function(data){
                                 'set the data'
                                 if(.validateColumns(colnames(data))){
                                   transitiondata <<- data[,.formalAllDef]
                                   precdata <<- unique(transitiondata[,.formalPrecDef])
                                   
                                 }
                               },
                               
                               getTransIntensity = function() {
                                 'matrix with transitions intensities'
                                 transIntensities <- dcast(transitiondata, ModifiedSequence + z + IonType + FragZ + Sequence ~ Filename , value.var="Intensity")
                                 return(transIntensities)
                               },
                               
                               getPrecursorIntensity=function(FUN=sum){
                                 transIntensities <- getTransIntensity()
                                 key <- paste(transIntensities$ModifiedSequence, ModifiedSequence$z, sep="_")
                                 res <- by(transIntensities,INDICES=key,FUN=FUN)
                                 do.call("rbind",res)
                               },
                               
                               getRT = function() {
                                 'matrix with precursor retention times'
                                 transRT = dcast(precdata, ModifiedSequence + z + Sequence  ~ Filename , value.var="RT")
                               },
                               
                               getScore = function() {
                                 'matrix with precursor scores'
                                 transScore = dcast(precdata, ModifiedSequence + z + Sequence ~ Filename , value.var="Score")
                                 return(transScore)
                               },
                               getMZ = function() {
                                 'matrix with precursor mz'
                                 transMz = dcast(precdata, ModifiedSequence + z + Sequence  ~ Filename , value.var="mz")
                                 return(transMz)
                               },
                               getDecoy= function(){
                                 'matrix with decoy information'
                                 transDecoy = dcast(precdata, ModifiedSequence + z + Sequence  ~ Filename , value.var="Decoy")
                                 return(transDecoy)
                               },
                               getFDR = function(){
                                 'compute FDR for dataset'
                                tmp <-data.frame(Protein = precdata$Protein,Decoy = precdata$Decoy,Score = precdata$Score)
                                tmp <- tmp[tmp$Score < 2,] # not sure how to treat requant values in this context
                                fdr <- (sum(tmp$Decoy) / length(tmp$Protein))
                                return(fdr)
                               },
                               plotFDR = function(log=""){
                                 'plot FDR versus Score'
                                 tmp <-data.frame(Protein = precdata$Protein,Decoy = precdata$Decoy,Score = precdata$Score)
                                 tmp <- tmp[tmp$Score < 2,] # not sure how to treat requant values in this context
                                 tmp <- tmp[order(tmp$Score),]
                                 plot(tmp$Score,cumsum(tmp$Decoy) / nrow(tmp) * 100 ,type="l",xlab="Score", ylab="FDR",log=log)
                               }
                               
                             ))


