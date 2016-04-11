MsExperiment <- setRefClass("msExperiment",
                            fields = list(
                              data = "data.frame",
                              #peptideIntensity = "data.frame",
                              #fragmentIntensity = "data.frame",
                              #mScore = "data.frame",
                              #rt = "data.frame",
                              #mz = "data.frame"

                            ),
                            methods = list(
                              setData <- function(data){data = data}
                            )
)


msexp <- MsExperiment$new(data=NULL)
msexp$data


class(msexp)
attributes(msexp)
msexp$setData("xx")
#person <- setRefClass("msExperiment")

person$methods(setData = function(data){data = data})
person$setData("xx")
