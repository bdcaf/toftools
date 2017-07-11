library(baseline)
library(zoo)

source('R/tofReader.R')
source('R/tofSparse.R')

#tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
tofA.h5 <- 'testdata/2017.06.22-11h13m23s ca valve open.h5'
tofB.h5 <- 'testdata/2017.06.20-15h04m22s Aurelio valve off.h5'
tofA <- tofH5(tofA.h5)
tofB <- tofH5(tofB.h5)
sumSpecA <- smooth(sumSpec.TofH5(tofA))
sumSpecB <- smooth(sumSpec.TofH5(tofB))

prep_spec <- function(sspec, atof, wid=30){
  sats <- find_saturated(sspec, N=atof$indexhelp$N)
  s1 <- dense_remove_sat(sspec, sats)
}

clA <- prep_spec(sumSpecA, tofA)
clB <- prep_spec(sumSpecB, tofB)

orientSpec <- smooth(log1p(clA))
toOrient <- smooth(log1p(clB))
plot(seq_along(orientSpec), orientSpec, type='l', xlim=c(11.1e4,13.5e4))
lines(toOrient, col='red')

#assure_range <- function(f)
   #function(x) pmax.int(pmin.int(f(x),length(normSpec)),1)

optim_generator <- function(orientSpec, wf){
  reference <- orientSpec/sqrt(sum(orientSpec^2))
  function(a, dense_spec){
	warp_fun <- wf(a)
	warped <- warp_dense( dense_spec, warp_fun)
	#warped$v <- log1p(warped$v)
	-cor.full.full(warped, reference) 
  }
}


#refSpec <- baseline(matrix(orientSpec, nrow=1), method= 'rollingBall', wm=100, ws=100)
refSpec <- baseline(matrix(orientSpec, nrow=1), method= 'rollingBall', wm=100, ws=100)@corrected[1,] # gut und generell
checkSpec <- baseline(matrix(toOrient, nrow=1), method= 'rollingBall', wm=100, ws=100)@corrected[1,] # gut und generell

refSpec <- pmax(refSpec,0)
checkSpec <- pmax(checkSpec,0)
wid <- 100
r2 <- stats::filter(as.numeric(refSpec), rep(1/wid, wid))

plot(refSpec, type='l', xlim=c(12.1e4,12.8e4))
lines(r2, col='green')
lines(checkSpec, col='red')

optim_fun <- optim_generator(r2, warp0)

warp0 <- function(a) 
  function(x) a[[1]]+x # + a[[2]]*x # + a[[3]]*x^2
opt_res <- optim( 0, optim_fun, gr=NULL, checkSpec, hessian = F, method='Brent', lower=-1000, upper=1000)
warped <- warp_dense( checkSpec, warp0(opt_res$par))
plot(seq_along(refSpec), refSpec, type='l', xlim=c(12.1e4,13.5e4))
lines(seq_along(checkSpec), (checkSpec), col='red')
with(warped, lines(starts:ends, v, col='blue'))
#Rprof()
#summaryRprof('work/profile')
warp1 <- function(a) 
  function(x) a[[1]] + a[[2]]*x # + a[[3]]*x^2
optim_fun1 <- optim_generator(refSpec, warp1)
#opt_res1 <- optim( c(10,1), optim_fun1, gr=NULL, checkSpec, hessian = F)
opt_res1 <- optim( c(0,1), optim_fun1, gr=NULL, checkSpec, hessian = F, 
				  method = 'L-BFGS-B', 
				  lower = c(-1000, 0.7),
				  upper = c(1000, 1.3))
warped <- warp_dense( checkSpec, warp0(opt_res$par))
plot(seq_along(refSpec), refSpec, type='l', xlim=c(12.1e4,13.5e4))
lines(seq_along(checkSpec), checkSpec, col='red')
with(warped, lines(starts:ends, v, col='blue'))
