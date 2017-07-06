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
sparse_spec <- function(full.wave, lower=0, minlen=1, max_gap = 100){
  tab <- data.table(id=seq_along(full.wave), wave=full.wave)[wave>lower]
  tab[, gap:= (id - shift(id, fill=0))]
  tab[, new_block := gap>max_gap]
  tab[, id_block := cumsum(new_block)]
  tab[new_block==F, 
	  v := Vectorize(function(gap,wave) c(rep(0,gap-1L),wave))(gap,wave)
  	  ]
  tab[new_block==T, v:=as.list(wave)]

  fullT <- tab[,.(starts=min(id), 
  				  ends=max(id), 
  				  v=list(Reduce(c,v))), 
  			   by=id_block]
  fullT[, len := ends-starts+1]
  fullT[ends>starts]
  return( fullT[len>=minlen])
}

simplify_sparse <- function(spec_pre, max_gap=10L){
  spec_pre[, gap := starts - shift(ends, fill = 0)][, inds := cumsum(gap > max_gap)]
  # close gaps but only within groups
  DT0 <- spec_pre[between(gap, 2L, max_gap), 
				  .(starts = starts - (gap - 1L), 
					ends = starts - 1L, 
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
  lastcu <- last(cuv)
  vc <- c(0,cuv,lastcu)
  i_new <- seq(from = floor(i_trans[[1]]), 
  			   to = ceiling(last(i_trans)))
  v_new <- diff( c(0,approx(i_trans, vc, 
  							i_new, yleft=0, yright=lastcu)$y))
  list(starts=i_new[1],  ends=i_new[length(i_new)], v=as.vector(v_new))
}
warp_spec <- function(semi_sparse, warpfun){
  vec_warp <- Vectorize( function(s,e,v) warp_line(warp_fun,s,e,v), SIMPLIFY=F)

  v2 <- function(s,e,v) mapply( function(x,y,z) warp_line(warp_fun,x,y,z),s,e,v)
  reorderer <- function(x,i) sapply(tmp, function(x) x[[i]])
  order_vec <- function(tmp) lapply(1:3, function(a) reorderer(tmp,a))
  vec_order_warp <- function(s,e,v) order_vec(vec_warp(s,e,v))
  warped <- semi_sparse[,vec_order_warp(starts,ends,v)]
  colnames(warped) <- c('starts','ends','v')

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
spsp <- sparse_spec(full.wave, lower=0, minlen=10, max_gap=30)
})
spsp[len>1000]
hist(spsp[len<100]$len)
system.time({
simplified <- simplify_sparse(spsp, max_gap=50L)
})

spt <- sparseTof(aSpec)
totalSpec <- sumSpec.TofH5(myTof) # sparsity makes no sense here
s

cor.semisparse.full <- function(){
  refSpec <- totalSpec
  refEnergy <- sum(totalSpec^2)
  spsp[, energy := Vectorize(function(x) sum(x^2))(v)]
  spsp[, cor := Vectorize( function(st,en,v) v %*% refSpec[st:en])(starts, ends,v)]
  agg <- spsp[, .(total_energy=sum(energy), total_sp=sum(cor)), by=NULL]
  with(agg, total_sp/sqrt(total_energy)/sqrt(refEnergy))
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
