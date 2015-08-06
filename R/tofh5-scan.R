#' Read TOF scan at specified index
#' @export
#' @docType methods
#' @rdname read.scan-methods
#' @exportMethod read.scan
setGeneric("read.scan", function(object, index) {
  standardGeneric("read.scan")
})


#' @rdname read.scan-methods
setMethod('read.scan', signature(object = 'TofH5', index='numeric'), function(object, index) {
  read.spec.ind(object@.tofblock, object@.indexhelp, index)
})
