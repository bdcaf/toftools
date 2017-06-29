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

  tmp <- 
  	data_frame(len = renc$lengths, 
  			   g0  = renc$values) %>%
	mutate( ends = cumsum(len),
		    starts = 1 + c(0, ends[-nrow(.)]) 
		    )

  small0 <- tmp %>% filter(g0==F, len<3)

  only_sig <- tmp %>% filter(g0) %>% select(starts, ends)
  condensed_sig <- TibAccRed(mergeFun, only_sig, max_dif=20)

  extract_ts <- function(starts, ends) list(full.wave[starts:ends])
  condensed_sig %>% rowwise() %>% mutate(v=extract_ts(starts,ends))
}

#' revert sparse spectrum to dense
#'
#' @description Turns a sparse spectrum into a dense vector.  Mostly
#' useful for quick testing.
#' @param sparse.spec, sparse spectrum
#' @param full.len, optional full length of dense spectrum
#' @return vector
#' @export
sparse2dense <- function(sparse.spec, full.len=398999){
  out <- vector(mode='numeric', length=full.len)
  sparse.spec %>% rowwise() %>%
  	do( x=with(., out[starts:ends] <<- c(unlist(v))))
  
  out
}
