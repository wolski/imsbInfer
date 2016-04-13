.formalAllDef <-c("z", "RT",  "ModifiedSequence",
                 "Sequence", "Filename", "Decoy",
                 "mz", "Score", "IonType",
                 "FragZ", "Intensity")

.formalPrecDef <- c("z", "RT",  "ModifiedSequence",
                   "Sequence", "Filename", "Decoy",
                   "mz", "Score")


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
#' huhu <- msTransitions()
#' huhu$columnsAll
#' huhu$columnsPrecursor
#' huhu$transitiondata
#' file = file.path(path.package("imsbInfer2"),"extdata/E1603291025_feature_alignment_requant.tsv.gz")
#' data <- read_tsv(file,col_names = TRUE)
#' data <- prepareOpenSwathData(data)
#' head(data)
#' huhu$setData(data)
#' intensity <- huhu$getTransIntensity()
#' dim(intensity)
#' rt <-huhu$getRT()
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
                                 transIntensities = dcast(transitiondata, ModifiedSequence + z + IonType + FragZ + Sequence ~ Filename , value.var="Intensity")
                                 return(transIntensities)
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
                               }
                               
                             ))


