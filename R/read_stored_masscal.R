#' reads the mass calibration from a file
#' @import rhdf5
#' @useDynLib rhdf5
#' @param file location of h5 file
#' @return list of conversion functions and parameters
read_stored_masscal <- function(file){
  fid1 <- H5Fopen(file)
  mcfun <- stored_mass_cal.tof_h5(fid1)
  H5Fclose(fid1)
  return(mcfun)
}
