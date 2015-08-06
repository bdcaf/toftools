#' open TOF file
#' 
#' required for reading tof data
#' @export
#' @docType methods
#' @rdname tofh5-methods
setGeneric("open", function(object) {
  standardGeneric("open")
})

#' @rdname tofh5-methods
#' @aliases tofh5,open,open-method
setMethod('open', signature(object = 'TofH5'), function(object) {
			object@.meas.reader <- H5Fopen(object@file.name)
			object@.tofblock <- get.raw.tofblock(object@meas.reader)
			object@.indexhelp <- tof.indexhelp(object@.tofblock)
})


#' close TOF file
#' 
#' usually not required
#' @export
#' @docType methods
#' @rdname tofh5-methods
setGeneric("close", function(object) {
  standardGeneric("close")
})

#' @rdname tofh5-methods
#' @aliases tofh5,close,close-method
setMethod('close', signature(object = 'TofH5'), function(object) {
			H5Fclose(object@meas.reader)
			object@.meas.reader <- NULL
			object@.tofblock <- NULL
			object@.indexhelp <- NULL
})


