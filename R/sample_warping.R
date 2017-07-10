source('R/tofReader.R')
source('R/tofSparse.R')
library(zoo)
#tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
tof.h5 <- 'testdata/2017.06.22-11h13m23s ca valve open.h5'
myTof <- tofH5(tof.h5)
aSpec <- readInd.TofH5(myTof,10)
full.wave <- aSpec
system.time({
spsp <- semisparse_spec(full.wave, lower=0, minlen=10, max_gap=30)
})
#system.time({
#simplified <- simplify_sparse(spsp, max_gap=50L)
#})

warp0 <- function(a) 
  function(x) a[[1]] + a[[2]]*x + a[[3]]*x^2

totalSpec <- sumSpec.TofH5(myTof)
refSpec <- smooth(totalSpec, twiceit=T)
#refSpec <- rollmean(totalSpec, 1000)
#refSpec <- totalSpec
#refSpec <- totalSpec
refEn <- sum(refSpec^2)

normSpec <- refSpec/sqrt(refEn)

assure_range <- function(f)
   function(x) pmax(1,pmin(f(x),length(normSpec)))



optim_fun <- function(a,spspec){
  wf <- assure_range(warp0(a))
  warped <- warp_spec(spspec, wf)

  #normfac <- sum(a^2) * 0.00001
  out <- -cor.semisparse.full(warped, normSpec) 
  #print(out)

  out
}

i<-1
a <- c(-10,1, -1.2/length(normSpec))
cspec <- readInd.TofH5(myTof, i)
spspec <- semisparse_spec(cspec, lower=0, minlen=10, max_gap=30)
startV <- c(0,1,0)
#Rprof('work/profile')
opt_res <- optim( startV, optim_fun, gr=NULL, spspec, 
				 method='BFGS',
				 #method='SANN',
				 #method='L-BFGS-B',
				 #lower = c(-Inf,0,-1.2/length(normSpec)),
				 #upper = c(Inf, 1.2, 1.2/length(normSpec)),
				 hessian = F)
#Rprof()
#summaryRprof('work/profile')

warp_par <- function(i){
  cspec <- readInd.TofH5(myTof, i)
  spspec <- semisparse_spec(cspec, lower=0, minlen=10, max_gap=30)
  opt_res <- optim(c(0,1,0), optim_fun, NULL, spspec, method='BFGS')
  opt_res$par
}

i_sel <- floor(seq(from =1 , to=myTof$indexhelp$N, length.out=20))
system.time( ww <- lapply(i_sel, warp_par))
offsets <- sapply(ww, function(x) x[[1]])
slopes <- sapply(ww, function(x) x[[2]])
plot(offsets)
plot(slopes)
