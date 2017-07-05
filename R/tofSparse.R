library(dplyr)
library(data.table)



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
sparse_spec <- function(full.wave, lower=0, minlen=1){
  vex <- Vectorize(function(st,en) (full.wave[st:en]), SIMPLIFY=F)

  idx <- full.wave > lower
  renc <- rle(as.vector(idx))

  edt <- data.table(lengths = renc$lengths, g0 = renc$values)
  edt[,ends := cumsum(lengths)][, starts := shift(ends, fill=0)+1]

  e2 <- edt[g0 == T & lengths >= minlen ,,]
  setkey(e2, starts, ends)
  e2[, c('g0','lengths'):=NULL]
  e2[, v := (vex(st=starts, en=ends))]
  e2
}

simplify_sparse <- function(spec_pre, max_gap=10L){
  spec_pre[, gap := starts - shift(ends, fill = 0)][, inds := cumsum(gap > max_gap)]
  # close gaps but only within groups
  DT0 <- spec_pre[between(gap, 2L, max_gap), .(starts = starts - (gap - 1L), ends = starts - 1L, 
										 v = Vectorize(rep.int)(0L, gap - 1L), gap, inds)]
  # bind rowwise (union in SQL), setkey on result to maintain sort order, 
  # remove column gap as no longer needed
  DT2 <- setkey(rbind(spec_pre, DT0), starts, ends)[, gap := NULL][]
  # aggregate groupwise, pick min/max, combine lists
  result <- DT2[, .(starts = min(starts), ends = max(ends), v = list(Reduce(c, v))), by = inds]
  # alternative code: pick first/last
  result <- DT2[, .(starts = first(starts), ends = last(ends), v = list(Reduce(c, v))), by = inds]
  result
}

warp_fun <- function(x) 0.99*x + 0.0000001*x^2 - 0.3
starts <- 75872
ends <- 75881
v <-c(1.00053060054779, 2.00164222717285, 2.00208926200867, 5.00980424880981,
4.00963306427002, 5.01561975479126, 0, 4.01455307006836, 3.01183128356934,
1.00410747528076)
semi_sparse <- simplified[1:5,]
warp_line <- function(warp_fun, starts, ends, v){
  i_trans <- warp_fun((starts-1):(ends+1))
  cuv <- cumsum(v)
  lastcu <- cuv[length(cuv)]
  vc <- c(0,cuv,lastcu)
  i_new <- seq(from = floor(i_trans[[1]]), 
  			   to = ceiling(i_trans[length(i_trans)]))
  v_new <- diff( c(0,approx(i_trans, vc, i_new, yleft=0, yright=lastcu)$y))
  list(i_new[1],  i_new[length(i_new)], v_new)
}
warp_spec <- function(semi_sparse, warpfun){
  vec_warp <- Vectorize( function(s,e,v) warp_line(warp_fun,s,e,v), SIMPLIFY=F)
  vec_warp(c(1:3), c(4:6), list(c(1:4), c(2:5), c(3:6)))
  wip <- semi_sparse[,c('an', 'nv', 'nn') := vec_warp(starts,ends,v), with=F ][]

}


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

tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
myTof <- tofH5(tof.h5)
aSpec <- readInd.TofH5(myTof,10)
full.wave <- aSpec
system.time({
spsp <- sparse_spec(full.wave, lower=0, minlen=3)
})
system.time({
simplified <- simplify_sparse(spsp, max_gap=50L)
})
spt <- sparseTof(aSpec)

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
