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

lss <- function(ss, col='black') lines(x=seq_along(ss), y=ss, col=col)

#plot.new()
# show shift
plot(seq_along(sumSpecA), sumSpecA, type='l', xlim=c(9.7e4,10.5e4))
#lss(sumSpecA)
lss(sumSpecB, col='red')


prep_spec <- function(sspec, atof, wid=30){
  sats <- find_saturated(sspec, N=atof$indexhelp$N)
  s1 <- dense_remove_sat(sspec, sats)
}

clA <- prep_spec(sumSpecA, tofA)
clB <- prep_spec(sumSpecB, tofB)

plot(seq_along(clA), clA, type='l', xlim=c(11.1e4,13.5e4))
lss(clB, col='red')


refSpec <- clA
refEn <- sum(refSpec^2)
normSpec <- refSpec/sqrt(refEn)

wid <- 30
orientSpec <- rollmean(log1p(clA),wid)
orientSpec <- orientSpec/sqrt(sum(orientSpec^2))
checkSpec <- log1p(clB)
checkSpec <- checkSpec/sqrt(sum(checkSpec^2))
plot(seq_along(orientSpec), orientSpec, type='l', xlim=c(11.1e4,13.5e4))
lss(checkSpec, col='red')

assure_range <- function(f)
   function(x) pmax.int(pmin.int(f(x),length(normSpec)),1)

optim_generator <- function(orientSpec, wf)
  function(a, dense_spec){
	warp_fun <- assure_range(wf(a))
	warped <- warp_dense( dense_spec, warp_fun)
	warped$v <- log1p(warped$v)
	-cor.full.full(warped, orientSpec) 
  }

warp0 <- function(a) 
  function(x) a[[1]]+x # + a[[2]]*x # + a[[3]]*x^2

optim_fun <- optim_generator(orientSpec, warp0)

opt_res <- optim( 0, optim_fun, gr=NULL, clB, hessian = F, method='Brent', lower=-1000, upper=1000)
warped <- warp_dense( clB, warp0(opt_res$par))
warped$logv <-log1p(warped$v)
plot(seq_along(orientSpec), orientSpec, type='l', xlim=c(11.1e4,13.5e4))
lines(seq_along(clB), log1p(clB)/sqrt(sum(log1p(clB)^2)), col='red')
with(warped, lines(starts:ends, logv/sqrt(sum(logv^2)), col='blue'))
#Rprof()
#summaryRprof('work/profile')
warp1 <- function(a) 
  function(x) a[[1]] + a[[2]]*x # + a[[3]]*x^2
optim_fun1 <- optim_generator(orientSpec, warp1)
opt_res1 <- optim( c(0,1), optim_fun1, gr=NULL, clB, hessian = F, 
				  method = 'L-BFGS-B', 
				  lower = c(-10000, 0.5),
				  upper = c(10000, 1.5))

