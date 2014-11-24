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