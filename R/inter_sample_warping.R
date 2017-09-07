library(BB)
library(baseline)
library(zoo)

source('R/tofReader.R')
source('R/tofSparse.R')

#tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
tofA.h5 <- 'testdata/2017.06.22-11h13m23s ca valve open.h5'
tofB.h5 <- 'testdata/2017.06.20-15h04m22s Aurelio valve off.h5'
tofA <- tof_h5(tofA.h5)
tofB <- tof_h5(tofB.h5)
sum_specA <- smooth(sum_spec.tof_h5(tofA))
sum_specB <- smooth(sum_spec.tof_h5(tofB))

prep_spec <- function(sspec, atof, wid=30){
  refSpec <- baseline(matrix(sspec, nrow=1), method= 'rollingBall', wm=100, ws=100)@corrected[1,]
  sats <- find_saturated(sspec, N=atof$indexhelp$N)
  s1 <- dense_remove_sat(sspec, sats)
}

clA <- prep_spec(sum_specA, tofA)
clB <- prep_spec(sum_specB, tofB)

orientSpec <- log1p(clA)
toOrient <- log1p(clB)
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
optim_generator_sp <- function(orientSpec, wf){
  reference <- orientSpec/sqrt(sum(orientSpec^2))
  function(a, sp_spec){
	warp_fun <- wf(a)
	warped <- warp_spec( sp_spec, warp_fun)
	#warped$v <- log1p(warped$v)
	-cor.semisparse.full(warped, reference) 
  }
}


#refSpec <- baseline(matrix(orientSpec, nrow=1), method= 'rollingBall', wm=100, ws=100)
refSpec <- baseline(matrix(orientSpec, nrow=1), method= 'rollingBall', wm=100, ws=100)@corrected[1,] # gut und generell
checkSpec <- baseline(matrix(toOrient, nrow=1), method= 'rollingBall', wm=100, ws=100)@corrected[1,] # gut und generell

refSpec[refSpec < 0.3] <- 0
checkSpec[checkSpec < 0.3] <- 0

sum(checkSpec>0)/length(checkSpec)

wid <- 99
r2 <- rollmean(refSpec, wid, align='center', fill=0)


plot(refSpec, type='l', xlim=c(12.1e4,12.4e4))
lines(r2, col='green')
lines(checkSpec, col='red')

warp0 <- function(a) 
  function(x) a[[1]]+x # + a[[2]]*x # + a[[3]]*x^2

#spcheck <- semisparse_spec(checkSpec)
optim_fun <- optim_generator(r2, warp0)
optim_fun_sp0 <- optim_generator_sp(r2, warp0)
optim_fun(10, checkSpec)
optim_fun_sp0(10, spcheck)

opt_res <- optim( 100, optim_fun, gr=NULL, checkSpec, 
				 control=list(trace=3, ndeps=10),
				 method='Brent', lower=-2000, upper=2000
				 )
opt_res
#opt_res <- optim( 100, optim_fun, gr=NULL, checkSpec, 
				 #method='SANN')
warped <- warp_dense( checkSpec, warp0(opt_res$par))
plot(seq_along(refSpec), refSpec, type='l', xlim=c(12.1e4,12.5e4))
lines(seq_along(checkSpec), (checkSpec), col='red')
with(warped, lines(starts:ends, v, col='blue'))
#Rprof()
#summaryRprof('work/profile')
#
warp1 <- function(a) 
  function(x) a[[1]] + a[[2]]*x # + a[[3]]*x^2
optim_fun1 <- optim_generator(r2, warp1)
opt_res1 <- optim( c(opt_res$par,1), optim_fun1, gr=NULL, checkSpec, hessian = F)
#opt_res1 <- optim( c(opt_res$par,1), optim_fun1, gr=NULL, checkSpec, 
				  #method = 'L-BFGS-B', 
				  #lower = c(-1000, 0.7),
				  #upper = c(1000, 1.3))
opt_res1
warped <- warp_dense( checkSpec, warp1(opt_res1$par))
plot(seq_along(refSpec), refSpec, type='l', xlim=c(12.1e4,12.5e4))
lines(seq_along(checkSpec), checkSpec, col='red')
with(warped, lines(starts:ends, v, col='blue'))

warp2 <- function(a) 
  function(x) a[[1]] + a[[2]]*x  + a[[3]]*x^2
optim_fun2 <- optim_generator(r2, warp2)
opt_res2 <- optim( c(opt_res1$par,0), optim_fun2, gr=NULL, checkSpec, hessian = F)
#opt_res2 <- optim( c(opt_res1$par,0), optim_fun2, gr=NULL, checkSpec, hessian = F, 
				  #method = 'L-BFGS-B', 
				  #lower = c(-1000, 0.7, 0),
				  #upper = c(1000, 1.3, 1e4))
warped <- warp_dense( checkSpec, warp2(opt_res2$par))
plot(seq_along(refSpec), refSpec, type='l', xlim=c(12.1e4,12.5e4))
lines(seq_along(checkSpec), checkSpec, col='red')
with(warped, lines(starts:ends, v, col='blue'))
