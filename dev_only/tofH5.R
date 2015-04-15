library(rhdf5)

#' @export
make.tof.h5 <- function(tof.name){
  obj <- list('file' = tof.name, .fid=c(), .isOpen=FALSE)
  class(obj) <- "TofH5"
  return(obj)
}

#' @export
open <- function(x) UseMethod("open")
#' @export
open.TofH5 <- function(obj) {
  obj$.fid <-H5Fopen(tof.h5)       
  obj$.isOpen=TRUE
  return(obj)
}

#' @export
close <- function(x) UseMethod("close")
#' @export
close.TofH5 <- function(obj) {
  H5Fclose(obj$.fid)  
  obj$.isOpen=FALSE
  return(obj)
}


read.waveforms <- function(obj){
  at <- H5Aopen(obj$.fid,"NbrWaveforms")
  waveforms <- H5Aread(at)
  H5Aclose(at)
  return(waveforms)
}


read.peak.data <- function(obj){
  gr <- H5Gopen(obj$.fid,"PeakData")
  pdh <- H5Dopen(gr,"PeakTable")
  peak.data <- H5Dread(pdh)
  H5Dclose(pdh)
  pdh <- H5Dopen(gr,"PeakData") # peak value in counts per extractions
  peak.value <- H5Dread(pdh)
  H5Dclose(pdh)
  return(peak.value)
}

# open <- function(x) UseMethod("open")
# 
# count <- function(x) UseMethod("count")
# add <- function(x) UseMethod("add")
# 
# 
# 
# add.TofH5 <- function(obj){
#   print(obj)
#   obj$num <- obj$num + 1
#   3
# }
