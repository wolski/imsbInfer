.conditionsColumns <- c("FileName", "Condition",  "Replicate", "Fraction", "Run" )

msAnnotation  <- 
  setRefClass("msAnnotation",
              fields = list( data = "list",
                             conditioncolumns="character",
                             experiment = function(x){
                               if(missing(x)){
                                 data$experiment
                               }
                               data$experiment <<- x
                             },
                             annotation = function(x){
                               print("in annotation")
                               if(missing(x)){
                                 query <- c("Select * from Annotation")
                                 return( dbGetQuery(.data$con,query) ) 
                               }
                             }
              ),
              methods = list(
                initialize = function(data=NULL,...) {
                  conditioncolumns <<- .conditionsColumns
                  data$experiment <<- data
                },
                setAnnotation = function(annotation){
                  if(!sum(conditioncolumns %in% colnames(annotation)) == length(conditioncolumns)){
                    warning("there is no column" , conditioncolumns[!(conditioncolumns %in% colnames(annotation))] )
                  }
                  filenames <- data$experiment$getFileName()
                  if(!sum(filenames %in% annotation$FileName) == length(filenames)){
                    warning("there is no annotation for file" , filenames[!(filenames %in% annotation$FileName)] )
                  }
                  tmp <- copy_to(data$experiment$.data, annotation, "Annotation", temporary = FALSE)
                  dbListTables(.data$con)
                })
  )
