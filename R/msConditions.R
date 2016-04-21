.conditionsColumns <- c("FileName", "Condition",  "Replicate", "Fraction", "Run" )

msAnnotation  <- 
  setRefClass("msAnnotation",
              fields = list( msExperiment = "msTransitionExperiment",
                             conditioncolumns="character",
                             
                             annotation = function(x){
                               print("in annotation")
                               if(missing(x)){
                                 query <- c("Select * from Annotation")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             }
              ),
              methods = list(
                initialize = function(msExperiment=NULL,...) {
                  conditioncolumns <<- .conditionsColumns
                  msExperiment <<- msExperiment
                },
                setAnnotation = function(annotation){
                  if(!sum(conditioncolumns %in% colnames(annotation)) == length(conditioncolumns)){
                    warning("there is no column" , conditioncolumns[!(conditioncolumns %in% colnames(annotation))] )
                  }
                  msExperiment$getFilenames()
                  tmp <- copy_to(msExperiment$.data, annotation, "Annotation", temporary = FALSE)
                  dbListTables(.data$con)
                })
  )
