library(rhdf5)
# library(dplyr)

#' reads integrated peaks as integrated by PTR-MS Viewer.
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
ions <- list(h3o = 21.0220875, no=29.99744, aceton=59.04914, isopren=69.06988, h5o2=39.03265, h7o3=55.03897)

#' calculate mass from TOF index
mass.calc <- function(pars, vec) {(( vec - pars[[1]] )/ pars[[2]] )^2 }

#' find location of ion maximum in TOF spectrum
find_ion_tof <- function (ion.mass = 21.0220875, preliminary.mass, curr.spec.line) {  
  ig <- preliminary.mass > ion.mass-0.3 & preliminary.mass < ion.mass+0.3  
  pos <- match(TRUE,ig) + which.max(spectrum[ig])
  ifelse(curr.spec.line[pos] > 30, pos, NA)
}

#' calibrate a single spec line
mass.calib.tof.single <- function (ions, prelim.mass, curr.spec.line){
  tofs <- sapply(ions, function(x) find_ion_tof(x, prelim.mass, curr.spec.line))
  out <- data.frame(mass = mass.calc(co,1:ncol(tmp)), signal=curr.spec.line)
}