# get individual scan
if (is.null(getGeneric("getScan"))){
setGeneric("getScan",
             function(object, index) standardGeneric("getScan"))
}

setMethod(f="getScan", 
          signature = signature("TofMeasurement","numeric"),
          definition=function(object, index=1){
            read.spec.ind(object@.tofBlock, object@.indexHelp,  index) 
          })
