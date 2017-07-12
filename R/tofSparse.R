library(data.table)

#' create a sparse representation of the full TOF scan
#'
#' @description For many functions it seems helpful to use a sparse
#' spectrum instead of the dense one stored. A usual breath TOF
#' spectrum of 400,000 indices has less than 1% values greater 0.
#' Use lower and minlen to control granularity.
#'
#' The sparse format is a simple data.table where each entry has three
#' fields: start the index of start point, rlen the length of signal
#' >0, and signal the actual values > 0.
#'
#' @param full.wave a full dense TOF spectrum
#' @param lower parameter setting the lower number of counts
#' @param minlen minimum lenght of series > 0 - 
#' @param max_gap to further simplify densly store gaps below this
#' number
#'
#' @return list of sparse regions where the spectrum is > `lower`	
#' @export
semisparse_spec <- function(full.wave, lower=0L, minlen=1L, max_gap = 30L){
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

#' helper function to simplify sparse spectra 
#'
#' @description simplifies spectra by closing gaps smaller than max_gap.
#' In the end it will not be usefull for correlation finding.  But it
#' was the function that led to a number of improvements.
#' @param spec_pre semisparse spectrum to simplify.
#' @param max_gap maximum gap width allowed
#' @return a potentially simplified spectrum
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

#' helper function for warping a single semisparse spectrum line
warp_line <- function(warp_fun, starts, ends, v, check=T){
  i_trans <- warp_fun((starts-1):(ends+1))
  if (check) {
	ditrans <- diff(i_trans)
	if (any(ditrans <= 0 ))  return(list(starts=NA,ends=NA, v=sum(v)))
  }
  cuv <- cumsum(v)
  lastcu <- last(cuv)
  vc <- c(0,cuv,lastcu)
  i_new <- seq(from = ceiling(i_trans[[1]]), 
  			   to = floor(last(i_trans)))
  v_new <- diff( c(0,approx(i_trans, vc, 
  							i_new, yleft=0, yright=lastcu)$y))
  list(starts=i_new[1],  ends=i_new[length(i_new)], v=as.vector(v_new))
}

#' helper to rearrange output of vectorize for data.table
seq_list <- function(x,i) sapply(x, function(x) x[[i]])
#' helper specifically sequencing warpline
seq_wf <- function(tmp) lapply(1:3, function(a) seq_list(tmp,a))

#' warps a complete semisparse spectrum by warpfun
#'
#' @description Does the lifting by recalculating a warped semisparse
#' spectrum.
#' @param semi_sparse semisparse spectrum
#' @warp_fun a function that recalculates the indices.  Should be
#' monotonic increasing
#' @return warped semisparse spectrum
warp_spec <- function(semi_sparse, warp_fun){
  vec_warp <- Vectorize( function(s,e,v) warp_line(warp_fun,s,e,v), SIMPLIFY=F)
  vec_order_warp <- function(s,e,v) seq_wf(vec_warp(s,e,v))
  warped <- semi_sparse[,vec_order_warp(starts,ends,v)]
  colnames(warped) <- c('starts','ends','v')
  warped
}

warp_dense <- function(dense_spec, warp_fun){
  warp_line(warp_fun, 1, length(dense_spec), dense_spec)
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

extr_ref <- function(ref,st,en){
  if (is.na(st) || is.na(en)) return(0)
  ran <- st:en
  ran[ran<1] <- NA
  ran[ran>length(ref)] <- NA
  res <- ref[ran]
  res[is.na(res)] <- 0
  res
}

spfun <- function(st,en,v, ref){
  sum(v * ref[st:en], na.rm=T)
}

#' calculates distance between spectra
#'
#' @description calculates cosine distance between a semisparse spectrum
#' and a reference spectrum.  Practically the reference is a dense
#' spectrum.
#' @param aSpec semisparse spectrum
#' @param refSpec reference spectrum
#' @param refEnergy energy of the reference spectrum may be
#' pre-calculated
#' @return cosine distance of spectra
#'
cor.semisparse.full <- function(aSpec, refSpec){
  aSpec[, energy := Vectorize(function(x) sum(x^2))(v)]
  aSpec[, sp := Vectorize( function(st,en,v) spfun(st,en,v, refSpec))(starts, ends,v)]
  agg <- aSpec[, .(total_energy=sum(energy), total_sp=sum(sp)), by=NULL]
  with(agg, total_sp/sqrt(total_energy))
}

cor.full.full <- function(fullWarp, refSpec){
  if (is.na(fullWarp$starts)) return(0)
  effective_range <- with(fullWarp, starts:ends)
  effective_range[effective_range<1] <- NA
  effective_range[effective_range>length(refSpec)] <- NA
  #effective_range[effective_range<1] <- NA
  #effective_range[effective_range>length(refSpec)] <- NA
  sp <- with(fullWarp, sum(refSpec[effective_range] * v, na.rm=T))
  total_energy <- with(fullWarp, sum(v^2,na.rm=T))
  sp[[1]]/sqrt(total_energy)
}

##' revert sparse spectrum to dense
##'
##' @description Turns a sparse spectrum into a dense vector.  Mostly
##' useful for quick testing.
##' @param sparse.spec, sparse spectrum
##' @param full.len, optional full length of dense spectrum
##' @return vector
##' @export
#sparse2dense <- function(sparse.spec, full.len=398999){
  #out <- vector(mode='numeric', length=full.len)
  #sparse.spec %>% rowwise() %>%
      #do( x=with(., out[starts:ends] <<- c(unlist(v))))
  
  #out
#}


find_saturated <- function(v, N=1, slope_crit= -300){
  vv <- v/N
  which(diff(vv)<slope_crit)
}

dense_remove_sat <- function(densev, sats, hide_range = -100:300){
  inds <- unique(sapply(sats, function(x) x+hide_range))
  densev[inds] = 0
  densev
}

sparse_remove_sat <- function(spspec, sats){
  sa1 <- sats[[1]]
  for (sa1 in sats)
	spspec <- spspec[starts>sa1 | ends<sa1] # = not (starts<sa1 & ends>sa1)
  spspec
}
