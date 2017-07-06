source('R/tofReader.R')
#tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
tof.h5 <- 'testdata/2017.06.22-11h13m23s ca valve open.h5'
myTof <- tofH5(tof.h5)
aSpec <- readInd.TofH5(myTof,10)
full.wave <- aSpec
system.time({
spsp <- sparse_spec(full.wave, lower=0, minlen=10, max_gap=30)
})
#system.time({
#simplified <- simplify_sparse(spsp, max_gap=50L)
#})

warp0 <- function(a) 
  function(x) a[[1]] + a[[2]]*x + a[[3]]*x^2

refSpec <- totalSpec
refEn <- sum(refSpec^2)

optim_fun <- function(a,cspec){
  wf <- warp0(a)
  warped <- warp_spec(cspec, wf)
  -cor.semisparse.full(warped, totalSpec, refEn)
}
warp_par <- function(i){
  cspec <- readInd.TofH5(myTof, i)
  spspec <- sparse_spec(cspec, lower=0, minlen=10, max_gap=30)
  opt_res <- optim(c(0,1,0), optim_fun, NULL, spspec)
  opt_res$par
}

i_sel <- floor(seq(from =1 , to=myTof$indexhelp$N, length.out=20))
system.time( ww <- lapply(i_sel, warp_par))
offsets <- sapply(ww, function(x) x[[1]])
slopes <- sapply(ww, function(x) x[[2]])
plot(offsets)
plot(slopes)
