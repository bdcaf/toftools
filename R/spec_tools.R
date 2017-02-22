#' create a sparse representation of the full TOF scan
#'
#' @description For many functions it seems helpful to use a sparse
#' spectrum instead of the dense one stored. A usual breath TOF
#' spectrum of 400,000 indices has less than 1% values greater 0.
#'
#' The sparse format is a simple list where each entry has three
#' fields: start the index of start point, rlen the length of signal
#' >0, and signal the actual values > 0.
#'
#' @param full.wave a full dense TOF spectrum
#' @param lower parameter setting the lower number of counts
#'
#' @return list of sparse regions where the spectrum is > `lower`	
#'
#' @export
sparse_spec <- function(full.wave, lower=0){
  idx <- full.wave > lower
  renc <- rle(as.vector(idx))
  startp <- cumsum( c(0, as.numeric(renc$lengths)))

  extract_ts <- function(starter, len)
	list(start = starter, 
		 rlen = len,
		 signal = full.wave[starter + seq_len(len)])

  dv <- data.frame(start=startp[seq_along(renc$lengths)], 
				   rlen = renc$lengths, 
				   greater0 = renc$values) %>%
				filter(rlen >3, greater0) 

  with(dv, mapply(extract_ts, start, rlen))
}
