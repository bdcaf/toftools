#' Calibrate TOF file mass axis
#' 
#' performs an elaborate mass calibration
#' @export
#' @docType methods
#' @rdname find.ion.block-methods
#' @exportMethod find.ion.block
setGeneric("find.ion.block", function(object,ions,...) {
  standardGeneric("find.ion.block")
})

#' @rdname find.ion.block-methods
setMethod('find.ion.block', signature(object = 'TofH5', ions = 'data.frame'), 
          function(object, ions, 
                   preliminary.coeff=init.mass.calib(object), n.samp=100, num.knots=3) {
            scl <- sample.scans(object, n.samp=n.samp)
            read.scan <- function(i) read.raw.scan(object,i)
            locate.ions.block_(ions,scl,read.scan,preliminary.coeff=preliminary.coeff)            
          })