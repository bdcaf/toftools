#' get just the individual scan (no mass axis included)
#' 
#' @description this uses a continous index instead of the buf-write notation
#' 
#' @param object the TofMeasurement object
#' @param index 1 based index of scan
#' 
#' @return MassSpectrum object from MALDIquant
#' @include TofMeasurementClass.R
#' @rdname tofmeasurement-methods
#' @name getScan
#' @export
if (is.null(getGeneric("getJustScan"))){
setGeneric("getJustScan",
             function(object, index) standardGeneric("getJustScan"))
}
#' @rdname tofmeasurement-methods
#' @export
setMethod(f = "getJustScan", 
          signature = signature("TofMeasurement","numeric"),
          definition = function(object, index=1){
            spec <- read_spec_at(object@.tofBlock, object@.indexHelp,  index)
          
            return(spec)
          })
