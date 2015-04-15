library(rhdf5)

setClass("TofH5", 
         representation(file.name='character', num='numeric', .Data='list'), 
         prototype(file.name='', num=1))
 
setGeneric("open", function(object) {
  standardGeneric("open")
})
setGeneric("add1", function(object) {
  standardGeneric("add1")
})

setMethod('open', signature(object = 'TofH5'), function(object) {
  print(object@file.name)
})


setMethod("add1", signature(object = "TofH5"), function(object) {
  object@num <- object@num +1
  return(object)
})
