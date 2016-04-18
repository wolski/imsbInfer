.conditionsColumns <- c("FileName", "Condition",  "Replicate", "Fraction", "Run" )

msConditions  <- setRefClass("msConditions",
                                             fields = list( conditions = "data.frame", conditioncolumns="character"),
                                             methods = list(
                                               initialize = function(...) {
                                                 require(reshape2)
                                                 conditioncolumns <<- .conditionsColumns
                                               }))
                                               