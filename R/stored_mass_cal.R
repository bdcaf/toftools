#' reads mass calibration from H5 file
#' @note This mass calibration is obtained bz TOF DAQ software. 
#' It uses the averaged spectrum to obtain the calibration.  There is a second 
#' mass calibration table in the file which uses averages of chunks.  The other one is just not used in this package.
#' @note This function is currently only used for convenience and uses the 
#' simples formula, hardcoded.   I am not sure if I will reuse it or even
#' bother to adjust for all the available formulas.
#' @param fid a h5f file handle
#' @return list of functions to convert between mass and index
#' @useDynLib rhdf5
stored_mass_cal.tof_h5 <- function(fid){
  gr <- h5readAttributes(fid, "FullSpectra")
  with(gr,
       list(
  to_mass = Vectorize(function(i) (
   (i - `MassCalibration p2`) / `MassCalibration p1`) ^ 2),
  to_index = Vectorize(function(m)
   `MassCalibration p2` + `MassCalibration p1` * sqrt(m)),
            coefs = c(p1 = `MassCalibration p1`,
                      p2 = `MassCalibration p2`)
             ))
}
