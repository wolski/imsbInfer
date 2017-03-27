# Hello, world!
#
# This is an example function named 'hello' 
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
library(reshape2)

topMedianFunc <- function(data, top = 3){
  top <- min(top, nrow(data))
  mat <- data[,3:ncol(data)]
  aggRowMean <- apply(mat, 1, median, na.rm = TRUE)
  mat <- mat[order( aggRowMean, decreasing = TRUE )[1:top], , drop=FALSE]
  c(nrPeptides = nrow(data), dimmat1 = dim(mat)[1], dimmat2 = dim(mat)[2],apply(mat,2,mean,na.rm=TRUE))
}

#intProtMat <- ddply(intMat, .(PG.ProteinAccessions), topFunc)
