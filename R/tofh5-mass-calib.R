#' Calibrate TOF file mass axis
#' 
#' is required before working with data
#' @export
#' @docType methods
#' @rdname mass.calib-methods
#' @exportMethod mass.calib
setGeneric("mass.calib", function(object) {
  standardGeneric("mass.calib")
})


#' @rdname mass.calib-methods
setMethod('mass.calib', signature(object = 'TofH5', ions = 'data.frame'), function(object, ions, n.samp=100) {
  scl <- sample.scans(object, n.samp=n.samp)
  locate.ions.block_(ions,scl,read.scan,preliminary.coeff)
  
})


#' Sample a number of scans equally over the measurments
#' 
#' @export
#' @docType methods
#' @rdname sample.scans-methods
#' @exportMethod sample.scans
setGeneric("sample.scans", function(object, ...) {
  standardGeneric("sample.scans")
})

#' @rdname sample.scans-methods
setMethod('sample.scans', signature(object = 'TofH5'), function(object, n.samp=100) {
  with(object@.indexhelp, 
       if(N > 2*n.samp) { floor(seq(from=1,to=N, length.out=n.samp))}
       else {  seq(from=1,to=N)}
  )
})