#' Read TOF scan at specified index
#' @export
#' @docType methods
#' @rdname read.raw.scan-methods
#' @exportMethod read.raw.scan
setGeneric("read.raw.scan", function(object, index) {
  standardGeneric("read.raw.scan")
})


#' @rdname read.raw.scan-methods
setMethod('read.raw.scan', signature(object = 'TofH5', index='numeric'), function(object, index) {
  read.spec.ind(object@.tofblock, object@.indexhelp, index)
})


#' Read TOF scan binned at provided masses
#' @export
#' @docType methods
#' @rdname read.bins-methods
#' @exportMethod read.bins
setGeneric("read.bins", function(object, ...) {
  standardGeneric("read.bins")
})


#' @rdname read.bins-methods
#' @useDynLib tofTools
setMethod('read.bins', signature(object = 'TofH5'), 
          function(object, mass.calibration, index, masses) {
  targets <- expand.grid(scan=index, ion=masses)  
  targets <- enrich.ions.block(targets)
  targets$idx <- predict(mass.calibration, newdata=targets)
  
  one.scan <- function(ct){
    sc <- read.raw.scan(object,ct[[1,'scan']])
    bin.vals<-read_bin_c(sc, ct$idx)
    
    with(ct,data.frame(scan=scan, ion=ion, bin=bin.vals))
  }
  
  group_by(targets, scan) %>% do(one.scan(.))
})