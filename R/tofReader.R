library(rhdf5)
library(dplyr)
library(tidyr)

#' reads integrated peaks as integrated by PTR-MS TOD-DAQ-Viewer.
#' 
#' WARNING: these must be first integrated in the software otherwise and error
#' will be thrown.
#' @export
#' @param tof.h5 filename of hdf5 file
#' @return dataframe with structure: rows = scans, cols = counts + information
#' @examples
#' tof.h5 <- file.path('data',"uncal.h5")
read.tof.peaks <- function( tof.h5 ){
  fid <-H5Fopen(tof.h5)
  at <- H5Aopen(fid,"NbrWaveforms")
  waveforms <- H5Aread(at)
  H5Aclose(at)
  pd1 <- H5Dopen(fid,"PeakData/PeakTable")
  peak.data <- H5Dread(pd1)
  pd2 <- H5Dopen(fid,"PeakData/PeakData") # peak value in counts per extractions
  peak.value <- H5Dread(pd2)
  
  peak.counts <- peak.value * waveforms[[1]]
  tmp <- apply(peak.counts, 1, as.vector)    
  # plot(tmp[1:200,40]) # plot H2O*H3O+
  peak.frame <- as.data.frame(tmp)
  colnames(peak.frame) <- peak.data$label
  peak.frame$file <- tof.h5
  return(peak.frame)
}

read.sum.spec <- function(fid){
  gr <- H5Dopen(fid, 'FullSpectra/SumSpectrum')
  H5Dread(gr)
}
read.mass.cals <- function(fid){
  gr <- H5Dopen(fid, 'FullSpectra/MassCalibration')
  H5Dread(gr)
}

#' approximate mass calibration
pars.approx <- list(intercept = 900,square_mass = 17600 )

#' sample ions for TOF calibration
ions <- c(h3o = 21.0220875, 
          no=29.99744, 
          o2 = 31.989281,
          aceton=59.04914,  
          h5o2=39.03265, 
          h7o3=55.03897,
          meth_formate = 61.028406,
          methacrolein = 71.04914)

#' calculate mass from TOF index
mass.calc <- function(pars, vec) {(( vec - pars[[1]] )/ pars[[2]] )^2 }

#' find location of ion maximum in TOF spectrum
find_ion_tof <- function (ion.mass = 21.0220875, preliminary.coeff, curr.spec.line, wid=0.3) {  
  mass.range <- round(preliminary.coeff$intercept + preliminary.coeff$square_mass*sqrt(ion.mass+ wid*c(-1,1)))
  ig <- mass.range[1]:mass.range[2]
  # plot(curr.spec.line[ig], type='l')
  pos <- ig[1]-1 + which.max(curr.spec.line[ig])
  ifelse(curr.spec.line[pos] > 10 & curr.spec.line[pos] < 1e6, pos, NA) # NA if either saturated or not enough
}

# #' calibrate a single spec line
# mass.calib.tof.single <- function (n.ions, preliminary.coeff, curr.spec.line) {
#   tofs <- sapply(n.ions, function(x) find_ion_tof(x, preliminary.coeff, curr.spec.line))
#   mf <- data.frame(sqmass=sqrt(n.ions), tofs = (tofs))
#   mod <- lm(tofs ~  sqmass, mf)
#   co <- coef(mod)
#   out <- data.frame(mass = mass.calc (co,1:length(curr.spec.line)), signal=curr.spec.line)
#   return(out)
# }

#'
mass.calib.coeff.single <- function (ions, preliminary.coeff, curr.spec.line) {
  tofs <- sapply(ions, function(x) find_ion_tof(x, preliminary.coeff, curr.spec.line))
  mf <- data.frame(sqmass=sqrt(ions), tofs = (tofs))
  mod <- lm(tofs ~  sqmass, mf)
  co <- coef(mod)
  out <- data.frame(intercept = co[[1]], square_mass = co[[2]])
  return(out)
}

#' calculates mass calibration coefficients for every scan in the flat tof array
#' @note this thkes about 30s for a measurement containing 3000 scans
#' @param tofblock from get.raw.tofblock
#' @param indexhelp from tof.indexhelp
#' @param preliminary.coeff some start values for coefficients
#' @param n.ions ion list to be used for calibration
#' @return data frame of calibration coefficients
#' @examples
#' tof.h5 <- file.path('data',"Ac just breath after C (2014-10-23T11h34m53_#).h5")
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' cr <- make.curr.tofreader(tofblock, indexhelp)
#' mc.table <- mass.calib.tof(tofblock, indexhelp)
#' @export
mass.calib.tof <- function(tofblock, indexhelp,
                           preliminary.coeff=list(intercept = 900,square_mass = 17650 ), 
                           n.ions = unlist(ions)){
  data.frame(scan = (1:indexhelp$N) ) %>% rowwise() %>% do( {
    curr.spec.line <- read.spec.ind(tofblock, indexhelp, .$scan)
    mc <- mass.calib.coeff.single(ions, preliminary.coeff=pars.approx, curr.spec.line)
    data.frame(scan=.$scan, intercept=mc[['intercept']], square_mass=mc[['square_mass']])
  } )
}

#' smoothes the mass calibration to avoid local jumps
#' @param mc.table raw mass calibration
#' @return mass calibration table after smoothing
#' @examples
#' mc.table <- mass.calib.tof(tofblock, indexhelp)
#' smooth.mc <- smooth.mass.cal(mc.table)
smooth.mass.cal <- function (mc.table) {
  mc.fit <- mc.table
  mc.fit$ind <- 1:nrow(mc.fit)  
  cc <- lm(cbind(intercept, square_mass) ~ ind, mc.fit)
  pv <- predict(cc,  data.frame(ind=1:nrow(mc.fit)))  
}

#' first step of reading - get the tof block
#' not export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' H5close(fid)
get.raw.tofblock <- function(fid) {
  gr <- H5Gopen(fid,"FullSpectra")
  H5Dopen(gr,"TofData")}

#' pre calculates some values
#' not export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
tof.indexhelp <- function(tofblock){
  h5spaceFile <- H5Dget_space(tofblock)
  dims <- H5Sget_simple_extent_dims(h5spaceFile)
  calc.indices <- expand.grid(buf=1:dims$size[[3]], write=1:dims$size[[4]])
  h5spaceMem <- H5Screate_simple(dims$size[[1]])
  list(N=nrow(calc.indices),
       h5spaceFile=h5spaceFile, 
       h5spaceMem=h5spaceMem, 
       calc.indices=calc.indices, 
       dims=dims$size)
}

#' reads the spec of a single tof scan
#' @export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' spec <- read.spec.ind(tofblock, indexhelp, 40)
#' plot(spec,type='l')
read.spec.ind <- function(tofblock, indexhelp, i){
  pos <- indexhelp$calc.indices[i,]
  H5Sselect_hyperslab(indexhelp$h5spaceFile, 
                      start = c(1,1,pos$buf,pos$write), 
                      count = c(indexhelp$dims[[1]],1,1,1) )                          
  H5Dread(h5dataset = tofblock, 
          h5spaceFile = indexhelp$h5spaceFile, 
          h5spaceMem = indexhelp$h5spaceMem )
}

#' helper to wrap reader call
#' not export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' cr <- make.curr.tofreader(tofblock, indexhelp)
#' spec5 <- cr(5)
make.curr.tofreader <- function(tofblock, indexhelp){
  function(i){read.spec.ind(tofblock,indexhelp,i)}
}


#' mass calibration function to simulate it as performed by breath view
#' @export
#' @param fid h5 file handle of the measurement
#' @param indexhelp shape of the data block
#' @return sorted data_frame where buf and write are from indexhelp,
#' while a and b are the coefficients for mass calibration 
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' mcv <- masscal_legacy_plain(fid, indexhelp)
masscal_legacy_plain <- function(fid, indexhelp){
  legacy_masscal <- read.mass.cals(fid)
  df <- as_data_frame(t(legacy_masscal)) %>%
	mutate(write = 1:nrow(.)) %>% rename(a=V1, b=V2)
  with(indexhelp, left_join(calc.indices, df, by=c('write')))
}

#' mass calibration function to improve it as performed by breath view
#' @export
#' @param fid h5 file handle of the measurement
#' @param indexhelp shape of the data block
#' @return sorted data_frame where buf and write are from indexhelp,
#' while a and b are the coefficients for mass calibration 
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' mcv <- masscal_legacy_plain(fid, indexhelp)
masscal_legacy_smooth <- function(fid, indexhelp){
  legacy_masscal <- read.mass.cals(fid)
  df <- as_data_frame(t(legacy_masscal)) %>% mutate(write = 1:nrow(.))
  a_smooth <- with(df, smooth.spline(write, V1))
  b_smooth <- with(df, smooth.spline(write, V2))
  ci <- indexhelp$calc.indices %>% 
	mutate(x = write - 0.5 + r*buf,
		   a = predict(a_smooth, x=x)$y,
		   b = predict(b_smooth, x=x)$y
		   )
  select(ci, buf, write, a,b)
}

#' make a mass axis based on resolution
#' no export
#' @param resolution of the mass axis in 1/r
#' @param mass.low starting point of mass axis (should be >0)
#' @param mass.high upper limit of mass axis
#' @return list of masses
make_massaxis <- function( resolution=1e4,
						   mass.low=10,
						   mass.high=200){
  log_mass_axis <- seq(from=log10(mass.low), 
					   to = log10(mass.high),
					   length.out=((mass.high-mass.low)*resolution)
					   )
  mass_axis <- 10^log_mass_axis
}

#' convert mass axis to indices in a specific spectrum
#' no export
#' @param calib containing a and b of mass calibration
#' @param target mass axis
#' @export vector of indices (non integer!)
massaxis2ind <- function(calib, mass_axis){
  with(calib,a*sqrt(mass_axis) + b)
}


#' extract resampled vector at mass axis
#' incoming counts are expected to arrive uniformly in measurement
#' interval
#' no export
#' @param spec vector of spectrum in raw format
#' @param calib current mass calibration
#' @mass_axis target mass axis
#' @return vector of resampled counts - not it is open to the left
#' with first value NA
spec2mav <- function(spec, calib, mass_axis){
  inds <- massaxis2ind(calib, mass_axis)
  cspec <- cumsum(spec)
  y2 <- approx(x=seq_along(cspec), y=cspec, xout=inds)$y
  mav <- c(NA, diff(y2))
}
