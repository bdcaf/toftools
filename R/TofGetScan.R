library(MALDIquant)

#' get individual scan
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
if (is.null(getGeneric("getScan"))){
setGeneric("getScan",
             function(object, index) standardGeneric("getScan"))
}
#' @rdname tofmeasurement-methods
#' @export
setMethod(f="getScan",
          signature = signature("TofMeasurement","numeric"),
          definition=function(object, index=1){
            spec <- read_spec_at(object@.tofBlock, object@.indexHelp,  index)
            a <- object@metaDataScan[['MassCalibration a']]
            b <- object@metaDataScan[['MassCalibration b']]
            mv <- .idx2mass(1:length(spec),a,b)

            s <- createMassSpectrum(mass=mv, intensity=as.numeric(spec),
                                    metaData=list(index = index))
            return(s)
          })
