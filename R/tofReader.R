library(rhdf5)
# library(dplyr)

#' reads integrated peaks as integrated by PTR-MS TOD-DAQ-Viewer.
#' 
#' WARNING: these must be first integrated in the software otherwise and error
#' will be thrown.
#' @export
#' @param tof.h5 filename of hdf5 file
#' @return dataframe with structure: rows = scans, cols = counts + information
read.tof.peaks <- function( tof.h5 ){
  fid <-H5Fopen(tof.h5)
  at <- H5Aopen(fid,"NbrWaveforms")
  waveforms <- H5Aread(at)
  H5Aclose(at)
  gr <- H5Gopen(fid,"PeakData")
  pdh <- H5Dopen(gr,"PeakTable")
  peak.data <- H5Dread(pdh)
  H5Dclose(pdh)
  pdh <- H5Dopen(gr,"PeakData") # peak value in counts per extractions
  peak.value <- H5Dread(pdh)
  H5Dclose(pdh)
  H5Fclose(fid)
  
  peak.counts <- peak.value * waveforms[[1]]
  tmp <- apply(peak.counts, 1, as.vector)    
  # plot(tmp[1:200,40]) # plot H2O*H3O+
  peak.frame <- as.data.frame(tmp)
  colnames(peak.frame) <- peak.data$label
  peak.frame$file <- tof.h5
  return(peak.frame)
}

#' approximate mass calibration
pars.approx <- list(intercept = 900,slope = 17600 )

#' sample ions for TOF calibration
ions <- list(h3o = 21.0220875, 
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
find_ion_tof <- function (ion.mass = 21.0220875, preliminary.coeff, curr.spec.line) {  
  ig <- round(preliminary.coeff$intercept + preliminary.coeff$square_mass*sqrt(ion.mass)) + c(-700:+700)  
  pos <- ig[1]-1 + which.max(curr.spec.line[ig])
  ifelse(curr.spec.line[pos] > 10 & curr.spec.line[pos] < 10000, pos, NA) # NA if either saturated or not enough
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
mass.calib.coeff.single <- function (n.ions, preliminary.coeff, curr.spec.line) {
  tofs <- sapply(n.ions, function(x) find_ion_tof(x, preliminary.coeff, curr.spec.line))
  mf <- data.frame(sqmass=sqrt(n.ions), tofs = (tofs))
  mod <- lm(tofs ~  sqmass, mf)
  co <- coef(mod)
  out <- data.frame(intercept = co[[1]], square_mass = co[[2]])
  return(out)
}

#' converts the raw multidimensional TOF array into flat structure
#' @param raw.tof.array multidimensional array as present in hdf5 file, as obtained by read.fullspectra.tof.h5
#' @export flat tof array
#' @examples
#'  flat.tof <- to.flat.tof(raw.tof)
to.flat.tof <- function(raw.tof.array) array(raw.tof.array, dim=c(dim(raw.tof.array)[1],prod(dim(raw.tof.array)[2:4])))

#' reads the full spectrum matrix from Tofwerks hdf5 file
#' @param tof.h5 filename of the hdf5 file
#' @export raw tof array (4 dimensional)
#' @examples
#'  raw.tof <- read.fullspectra.tof.h5(file.path('data',"uncal.h5"))
read.fullspectra.tof.h5 <- function (tof.h5) h5read(tof.h5,'FullSpectra/TofData')

#' reads the full spectrum matrix from Tofwerks hdf5 file and converts it to flat
#' @param tof.h5 filename of the hdf5 file
#' @export flat tof array (2 dimensional)
#' @examples
#'  flat.tof <- read.tof.fullspectra(file.path('data',"uncal.h5"))
#'
#' @export
read.tof.fullspectra <- function(tof.h5){
  raw.tof <- read.fullspectra.tof.h5(tof.h5)
  to.flat.tof(raw.tof)
}

#' calculates mass calibration coefficients for every scan in the flat tof array
#' @param flat.tof array as obtained from to.flat.tof
#' @param preliminary.coeff some start values for coefficients
#' @param n.ions ion list to be used for calibration
#' @return data frame of calibration coefficients
#' @examples
#'  mc.table <- mass.calib.tof(flat.tof)
#' @export
mass.calib.tof <- function(flat.tof, 
                           preliminary.coeff=list(intercept = 900,square_mass = 17650 ), 
                           n.ions = unlist(ions)){
  ap.a <- apply(flat.tof,2, function(x) mass.calib.coeff.single(n.ions = unlist(ions), 
                                                        preliminary.coeff=preliminary.coeff, 
                                                        curr.spec.line = x))
  df.a <- as.data.frame(matrix(unlist(ap.a), ncol=2, byrow=TRUE))
  colnames(df.a) <- c('intercept','square_mass')
  df.a
}

#' smoothes the mass calibration to avoid local jumps
#' @param mc.table raw mass calibration
#' @returns mass calibration table after smoothing
#' @examples
#' smooth.mc <- smooth.mass.cal(mc.table)
smooth.mass.cal <- function (mc.table) {
  mc.fit <- mc.table
  mc.fit$ind <- 1:nrow(mc.fit)  
  aa <- lm(intercept ~ ind, mc.fit)
  bb <- lm(square_mass  ~ ind , mc.fit) 
  
  res <- data.frame(
    ind = 1:nrow(mc.fit),
    intercept = predict(aa,  data.frame(ind=1:nrow(mc.fit))),
    square_mass  = predict(bb,  data.frame(ind=1:nrow(mc.fit)))
  )
}