library(dplyr)
library(purrr)

#' create a sparse representation of the full TOF scan
#'
#' @description For many functions it seems helpful to use a sparse
#' spectrum instead of the dense one stored. A usual breath TOF
#' spectrum of 400,000 indices has less than 1% values greater 0.
#' Use lower and minlen to control granularity.
#'
#' The sparse format is a simple list where each entry has three
#' fields: start the index of start point, rlen the length of signal
#' >0, and signal the actual values > 0.
#'
#' @param full.wave a full dense TOF spectrum
#' @param lower parameter setting the lower number of counts
#' @param minlen minimum lenght of series > 0 - 
#'
#' @return list of sparse regions where the spectrum is > `lower`	
#'
#' @export
sparse_spec <- function(full.wave, lower=0, minlen=0){
  idx <- full.wave > lower
  renc <- rle(as.vector(idx))
  startp <- cumsum( c(0, as.numeric(renc$lengths)))

  extract_ts <- function(starter, len)
	list(start = starter, 
		 rlen = len,
		 signal = full.wave[starter + seq_len(len)])

  calc_inds <- function(starter, len) starter + seq_len(length.out = len)
}

inpterp_sparse <- function(spa, func){
# TODO
}
corr_sparse <- function(spa, spb){
# TODO
}
