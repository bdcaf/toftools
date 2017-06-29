library(dplyr)
library(purrr)

ifun <- function(ll, inds)  ll[[1]] + ll[[2]]*inds + ll[[3]]*inds^2


one_interp <- function(lin, func){
  cums <- cumsum(c(0,0,lin[['v']],0,0))
  bins <- with(lin, (starts-2):(ends+2))
  nbins <- ifun(func, bins)
  tbins <- floor(nbins[[1]]):ceiling(nbins[[length(nbins)]])
  cout <- approx(nbins, cums, tbins, rule=2)
  out <- diff(cout$y)

  data_frame(starts2=tbins[[2]], ends2=tbins[[length(tbins)]], v2=list(out))
}

interp_sparse <- function(spa, func){
  spa %>% rowwise() %>% do(one_interp(.,func))
}

one_cor <- function(dense, lin){
  sum(with(lin, dense[starts2:ends2] * unlist(v2)))
}
corr_sparse <- function(dense, spa2){
  tmp <- spa2 %>% rowwise() %>% do(cor=one_cor(dense,.))
  vals <- unlist(tmp$cor)
  sum(vals, na.rm=T)
}

ofun <- function(dense, spa, func){
  spa2 <- interp_sparse(spa, func)
  corr_sparse(dense, spa2)
}

# even slower!!!
#concfun <- function(x) ofun(dense, spa, x)
#of <- optim(c(0,1,0), concfun )

