#' class to hold raw tof access
#' @export
#' @param tof.h5 file path to relevant H5 file
#' @return list of functions
#' @useDynLib tofTools
#' @examples
#' \dontrun{
#' tof.h5 <- file.path('data/',"2015.07.17-10h40m34 Ethanol deurated Karl .h5")
#' curr.h5 <- TofRaw(tof.h5)
#' }
TofRaw <- function(tof.h5){
  fid <-H5Fopen(tof.h5)
  tofblock <- get.raw.tofblock(fid)
  indexhelp <- tof.indexhelp(tofblock)
  mass.calib <- NULL
  cr <- make.curr.tofreader(tofblock, indexhelp)  
  
  read.scan <- function (i) read.spec.ind(tofblock, indexhelp, i)
  
  sampleScans <- function(n=100) {
    with(indexhelp, 
         if(N > 2*n) { floor(seq(from=1,to=N, length.out=n))}
         else {  seq(from=1,to=N)}
    )}
  
  mass.calibrate <- function(df.ion, n.samp=100,preliminary.coeff=list(intercept = 900,square_mass = 17650 ),
                            num.knots=3){
    scl <- sampleScans(n.samp)
    pos.block <- locate.ions.block_(df.ion,scl,read.scan,preliminary.coeff=preliminary.coeff)
    pos.block <- enrich.ions.block(pos.block)
    kn <- floor(seq(1,to=indexhelp$N, length.out=num.knots+2))[2:num.knots+1]  
    mass.calib <<- rlm(pos ~ ion + sq3.mass + sq5.mass + sq.mass * ns(scan,knots=kn), pos.block, method='MM')
  }
  
  structure(list(close = function() H5close(),
				 readIndex = cr, 
				 sampleScans = sampleScans,
				 num.scan = function() indexhelp$N,
				 read.scan = read.scan,
				 read.bin = function (i, bin.vec) {
				   sc <- read.spec.ind(tofblock, indexhelp, i)
				   readScales(sc, bin.vec)
				 },
				 mass.calibrate = mass.calibrate,
				 get.mass.calib = function() mass.calib),
			class='TofRaw'
			)  
}
