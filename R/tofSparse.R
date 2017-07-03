library(dplyr)
library(purrr)
library(S4Vectors)
library(XVector)
library(Matrix)


sparseTof <- function(full.wave){
  sp_spec <- sparse_spec(full.wave)
  this <- list(semi_sparse=sp_spec,
  			   length = length(full.wave))
  structure(this, class="sparseTof")
}

display = function(obj) {
  UseMethod("display", obj)
}

display.sparseTof = function(obj,...) {
  cat("A semi-sparse Tof data set.\n")
}

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
sparse_spec <- function(full.wave, lower=0, minlen=0){
  idx <- full.wave > lower
  renc <- rle(as.vector(idx))

  extract_ts <- function(starts, ends) list(full.wave[starts:ends])
  sparse_frame <- 
  	with(renc, data_frame(lens=lengths,g0=values)) %>%
  	mutate(ends=cumsum(lens), starts = lag(ends, default=0)+1) %>%
  	filter(g0) %>% select(starts, ends) %>% 
  	rowwise() %>% mutate(v=extract_ts(starts,ends), len=ends-starts+1) %>%
  	ungroup()
}

orig_frame <- sparse_frame %>% filter(len>1)
simp_frame <- simplify_semisparse(sparse_frame)

join_lines <- function(dfi) {
  if (nrow(dfi)==1) return(dfi)
  else 
  	with(dfi,{ 
  			 start <- starts[[1]]
			 end <- ends[[length(ends)]]
			 vals <- numeric(end-start+1)
			 add_val <- function(ddf)
			   with(ddf,{ 
					  vals[(starts-start+1) : (ends-start+1)] <<- 
						vals[(starts-start+1) : (ends-start+1)] + v })
			 dfi %>% rowwise() %>% do(tmp=add_val(.))
  			 data_frame(starts=start, ends=end, v=list(vals),
  			   len = end-start+1, inds = inds[[1]], join_pre=F)})
}

simplify_semisparse <- function(orig_frame, closeness = 10){
	orig_frame %>% 
	mutate( join_pre = lag(ends, default=0)+closeness >= (starts),
		    inds = cumsum(!join_pre)
		   )
	w2 <- wip %>% group_by(inds) %>% do(join_lines(.)) %>% ungroup()
}


tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
myTof <- tofH5(tof.h5)
aSpec <- readInd.TofH5(myTof,10)
spsp <- sparseSpec(aSpec)
sparseSpec <- function(dense_spec){
  semisparse <- sparse_spec(dense_spec)
  g0 <- sum(sapply(semisparse$v, length))
  structure(list(total_n=length(dense_spec),
  				 greater0 = g0,
  				 spec_table=semisparse), 
  			class="sparseSpec")
}

print.SparseSpec <- function(spsp){
  with(spsp, cat(paste0("Semi-sparse spec, ",total_n," entries, ",
  						greater0,"(",round(greater0/total_n*100,1),"%) >0\n")))
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
