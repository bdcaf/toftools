library(MALDIquant)

# get individual scan
if (is.null(getGeneric("getScan"))){
setGeneric("getScan",
             function(object, index) standardGeneric("getScan"))
}

setMethod(f="getScan", 
          signature = signature("TofMeasurement","numeric"),
          definition=function(object, index=1){
            spec <- read.spec.ind(object@.tofBlock, object@.indexHelp,  index)
            a <- object@metaDataScan[['MassCalibration a']]
            b <- object@metaDataScan[['MassCalibration b']]
            mv <- .idx2mass(1:length(spec),a,b)
            plot(mv,spec, type='l')
          })
