#' Read TOF scan at specified index
#' @export
#' @docType methods
#' @rdname read.raw.scan-methods
#' @exportMethod read.raw.scan
setGeneric("read.raw.scan", function(object, index) {
  standardGeneric("read.raw.scan")
})


#' @rdname read.raw.scan-methods
setMethod('read.raw.scan', signature(object = 'tof_h5', index='numeric'), function(object, index) {
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
setMethod('read.bins', signature(object = 'tof_h5'), 
          function(object, mass.calibration, index, masses) {
  targets <- expand.grid(scan=index, ion=masses)  
  targets <- enrich.ions.block(targets)
  targets$idx <- predict(mass.calibration, newdata=targets)
  read.prepared.bins(object, targets)
})


#' Read TOF scan at prepared targets 
#' @param targets list of targets with column idx specifying the relevant TOF indices and column index the relevant scan number
#' @export
#' @docType methods
#' @rdname read.prepared.bins-methods
#' @exportMethod read.prepared.bins
setGeneric("read.prepared.bins", function(object, ...) {
  standardGeneric("read.prepared.bins")
})

#' @rdname read.prepared.bins-methods
#' @examples
#' \dontrun{
#' tof.h5 <- file.path('testdata/Ac just exhale (2014-10-23T10h43m42_#).h5')
#' data(sample_ions)
#' ions <- sample_ions
#' 
#' object <- new("tof_h5", file.name=tof.h5)
#' object <- open(object)
#' 
#' scl <- sample.scans(object, n.samp=100)
#' 
#' ion.block <- find.ion.block(object,ions)
#' 
#' 
#' ion.block <- enrich.ib(ion.block)
#' head(ion.block)
#' num.knots <- 2
#' kn <- floor(seq(1,to=object@@.indexhelp$N, length.out=num.knots+2))[2:num.knots+1]
#' fit.mos <- rlm(pos ~ (isq + ion + sq) * ns(scan,knots=kn), ion.block, method='MM')
#' 
#' extr.data <- expand.grid(scan= 1:10, ion=seq(47,47.1,length.out=100))
#' extr.data <- enrich.ib(extr.data)
#' extr.data$idx <- predict(fit.mos, newdata = extr.data)
#' 
#' read.prepared.bins(object,targets=extr.data)
#' }
#' @useDynLib tofTools
setMethod('read.prepared.bins', signature(object = 'tof_h5'), 
          function(object, targets) {
            one.scan <- function(ct){
              sc <- read.raw.scan(object,ct[[1,'scan']])
              bin.vals<-read_bin_c(sc, ct$idx)              
              with(ct,data.frame(scan=scan, ion=ion, bin=bin.vals))
            }
            group_by(targets, scan) %>% do(one.scan(.))
          })
