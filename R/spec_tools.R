library(dplyr)
library(purrr)

canMerge <- function(sp1, sp2, max_dif = 10) { (sp2[[1,'starts']] - sp1[[1,'ends']]) < max_dif }
mergeSp <- function(sp1, sp2){ 
  sp1$ends <- sp2$ends
  sp1
}

mergeFun <- function(acc, cur, tib, max_dif){
  if (canMerge(acc, cur, max_dif)){
	tibout <- tib
	acc2 <- mergeSp(acc, cur)
  } else {
	tibout <- rbind(tib, acc)
	acc2 <- cur
  }
  list(acc = acc2, tib = tibout)
}

TibAccRed <- function(f, tibin, max_dif=10){
  init <- data_frame()
  acc <- tibin[1,]
  tibin <- tibin[-1,]
  
  while (nrow(tibin) > 0){
	cur <- tibin[1,]
	tibin <- tibin[-1,]

  	val <- f(acc, cur, init, max_dif)
	acc <- val$acc
	init <- val$tib
  }

  rbind(init, acc) 
}

#dd <- TibAccRed(mergeFun, only_sig)

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

sparse2dense <- function(sparse.list, full.len=398999){
  out <- vector(mode='numeric', length=full.len)
  sparse.list %>% rowwise() %>%
  	do( x=with(., out[starts:ends] <<- c(unlist(v))))
  
  out
}

inpterp_sparse <- function(spa, func){
# TODO
}
corr_sparse <- function(spa, spb){
# TODO
}



