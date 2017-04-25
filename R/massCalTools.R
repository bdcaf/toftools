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
  r <- 1/indexhelp$dims[[3]]
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
make_mass_axis <- function( resolution=1e4,
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
#' @param mass_axis mass axis
#' @export vector of indices (non integer!)
massaxis2ind <- function(calib, mass_axis){
  with(calib,a*sqrt(mass_axis) + b)
}

#' convert indices to local mass axis in a specific spectrum
#' no export
#' @param calib containing a and b of mass calibration
#' @param ind indices in spectra
#' @export vector of indices (non integer!)
ind2massaxis <- function(calib, ind){
  with(calib,((ind-b)/a)^2 )
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
  y2 <- approx(x=seq_along(cspec), 
  			   y=cspec, 
  			   xout=inds,
  			   method="linear")$y
  mav <- c(NA, diff(y2))
}
