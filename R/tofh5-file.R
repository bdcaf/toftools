#' open TOF file
#' 
#' open required for reading tof data
#' @export
#' @docType methods
#' @rdname open-methods
#' @exportMethod open
setGeneric("open", function(object) {
  standardGeneric("open")
})

#' @rdname open-methods
setMethod('open', signature(object = 'TofH5'), function(object) {
			object@.meas.reader <- H5Fopen(object@file.name)
			object@.tofblock <- get.raw.tofblock(object@.meas.reader)
			object@.indexhelp <- tof.indexhelp(object@.tofblock)
      object
})


#' close TOF file
#' 
#' close usually not required
#' @export
#' @docType methods
#' @rdname close-methods
#' @exportMethod close
setGeneric("close", function(object) {
  standardGeneric("close")
})

#' @rdname close-methods
setMethod('close', signature(object = 'TofH5'), function(object) {
			H5Fclose(object@.meas.reader)		
      
			object
})


