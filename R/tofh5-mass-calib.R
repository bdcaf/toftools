library(MASS)
library(splines)

#' Calibrate TOF file mass axis
#' 
#' performs an elaborate mass calibration
#' @export
#' @docType methods
#' @rdname mass.calib-methods
#' @exportMethod mass.calib
setGeneric("mass.calib", function(object,ions,...) {
  standardGeneric("mass.calib")
})

#' @rdname mass.calib-methods
setMethod('mass.calib', signature(object = 'TofH5', ions = 'data.frame'), 
          function(object, ions, 
                   preliminary.coeff=init.mass.calib(object), n.samp=100, num.knots=3) {
  scl <- sample.scans(object, n.samp=n.samp)
  read.scan <- function(i) read.raw.scan(object,i)
  ions.block <- locate.ions.block_(ions,scl,read.scan,preliminary.coeff=preliminary.coeff)
  ions.block <- enrich.ions.block(ions.block)
  kn <- floor(seq(1,to=object@.indexhelp$N, length.out=num.knots+2))[2:num.knots+1]  
  # object@.mass.calib <- 
  rlm(pos ~ ion + sq3.mass + sq5.mass + sq.mass * ns(scan,knots=kn), ions.block, method='MM')
#   object
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




#' Read inital guesses for mass calibration
#' 
#' read from the H5 file the calibrated value from Ionicon
#' @docType methods
#' @rdname init.mass.calib-methods
setGeneric("init.mass.calib", function(object,...) {
  standardGeneric("init.mass.calib")
})


#' @rdname init.mass.calib-methods
setMethod('init.mass.calib', signature(object = 'TofH5'), function(object) {
  gr <- H5Gopen(object@.meas.reader,"FullSpectra")
  ao1 <- H5Aopen(gr,'MassCalibration a')
  ao2 <- H5Aopen(gr,'MassCalibration b')
  a <- H5Aread(ao1)
  b <- H5Aread(ao2)
  
  H5Aclose(ao1)
  H5Aclose(ao2)
  list(square_mass=a, intercept=b)  
})
