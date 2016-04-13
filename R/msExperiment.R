formalAllDef <-c("z", "RT",  "ModifiedSequence",
                 "Sequence", "Filename", "Decoy",
                 "mz", "Score", "IonType",
                 "FragCharges", "Intensity")

formalPrecDef <- c("z", "RT",  "ModifiedSequence",
                   "Sequence", "Filename", "Decoy",
                   "mz", "Score")



## a simple editor for matrix objects.  Method  $edit() changes some
## range of values; method $undo() undoes the last edit.
msTransitions <- setRefClass("msTransitions",
                             fields = list( transitiondata = "data.frame",
                                            precdata = "data.frame",
                                            columnsAll="character",
                                            columnsPrecursor="character"),
                             methods = list(
                               initialize = function(...) {
                                 .formalAllDef <-c("z", "RT",  "ModifiedSequence",
                                                  "Sequence", "Filename", "Decoy",
                                                  "mz", "Score", "IonType",
                                                  "FragCharges", "Intensity")
                                 
                                 .formalPrecDef <- c("z", "RT",  "ModifiedSequence",
                                                    "Sequence", "Filename", "Decoy",
                                                    "mz", "Score")
                                 
                                 columnsAll <<- formalAllDef
                                 columnsPrecursor <<- formalPrecDef
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
                                 transRT = dcast(precData, ModifiedSequence + z + Sequence  ~ Filename , value.var="RT")
                                 
                               },
                               
                               getScore = function() {
                                 'matrix with precursor scores'
                                 transScore = dcast(precData, ModifiedSequence + z + Sequence ~ Filename , value.var="Score")
                                 return(transScore)
                               }
                               
                               getMZ = function() {
                                 'matrix with precursor mz'
                                 transMz = dcast(precData, ModifiedSequence + z + Sequence  ~ Filename , value.var="mz")
                                 return(transMz)
                               }
                               
                             ))

huhu <- msTransitions()
huhu$columnsAll
huhu$columnsPrecursor
huhu$rawdata

huhu$setData()
## add a method to save the object
mEdit$methods(
  save = function(file) {
    'Save the current object on the file
    in R external object format.
    '
    base::save(.self, file = file)
  }
)

tf <- tempfile()
xx$save(tf)



## Not run:
## Inheriting a reference class:  a matrix viewer
mv <- setRefClass("matrixViewer",
                  fields = c("viewerDevice", "viewerFile"),
                  contains = "mEdit",
                  methods = list( view = function() {
                    dd <- dev.cur(); dev.set(viewerDevice)
                    devAskNewPage(FALSE)
                    matplot(data, main = paste("After",length(edits),"edits"))
                    dev.set(dd)},
                    edit = # invoke previous method, then replot
                      function(i, j, value) {
                        callSuper(i, j, value)
                        view()
                      }))

## initialize and finalize methods
mv$methods( initialize =
              function(file = "./matrixView.pdf", ...) {
                viewerFile <<- file
                pdf(viewerFile)
                viewerDevice <<- dev.cur()
                dev.set(dev.prev())
                callSuper(...)
              },
            finalize = function() {
              dev.off(viewerDevice)
            })

## debugging an object: call browser() in method $edit()
xx$trace(edit, browser)

## debugging all objects from class mEdit in method $undo()
mEdit$trace(undo, browser)
## End(Not run)
