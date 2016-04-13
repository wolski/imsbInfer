## a simple editor for matrix objects.  Method  $edit() changes some
## range of values; method $undo() undoes the last edit.
mEdit <- setRefClass("mEdit",
                     fields = list( data = "matrix",
                                    edits = "list"),
                     methods = list(
                       edit = function(i, j, value) {
                         ## the following string documents the edit method
                         'Replaces the range [i, j] of the
                         object by value.
                         '
                         backup <-
                           list(i, j, data[i,j])
                         data[i,j] <<- value
                         edits <<- c(edits, list(backup))
                         invisible(value)
                       },
                       undo = function() {
                         'Undoes the last edit() operation
                         and update the edits field accordingly.
                         '
                         prev <- edits
                         if(length(prev)) prev <- prev[[length(prev)]]
                         else stop("No more edits to undo")
                         edit(prev[[1]], prev[[2]], prev[[3]])
                         ## trim the edits list
                         length(edits) <<- length(edits) - 2
                         invisible(prev)
                       },
                       show = function() {
                         'Method for automatically printing matrix editors'
                         cat("Reference matrix editor object of class",
                             classLabel(class(.self)), "\n")
                         cat("Data: \n")
                         methods::show(data)
                         cat("Undo list is of length", length(edits), "\n")
                       }
                     ))

xMat <- matrix(1:12,4,3)
xx <- mEdit()
xx$data

xx$edit(2, 2, 0)
xx
xx$undo()
mEdit$help("undo")
stopifnot(all.equal(xx$data, xMat))

utils::str(xx) # show fields and names of non-trivial methods

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
