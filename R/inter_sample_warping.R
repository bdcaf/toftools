library(magrittr)
library(rCharts)

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


prep_spec <- function(sspec, atof){
sats <- find_saturated(sspec, N=atof$indexhelp$N)
dense_remove_sat(sspec, sats)
}

clA <- prep_spec(sumSpecA, tofA)
clB <- prep_spec(sumSpecB, tofB)

plot(seq_along(clA), clA, type='l', xlim=c(11.1e4,13.5e4))
lss(clB, col='red')


warp0 <- function(a) 
  function(x) a[[1]]+x # + a[[2]]*x # + a[[3]]*x^2

refSpec <- clA
refEn <- sum(refSpec^2)
normSpec <- refSpec/sqrt(refEn)

assure_range <- function(f)
   function(x) pmax.int(pmin.int(f(x),length(normSpec)),1)

optim_generator <- function(wf)
  function(a, dense_spec){
	warp_fun <- assure_range(warp0(a))
	warped <- warp_dense( dense_spec, warp_fun)
	-cor.full.full(warped, normSpec) 
  }

warp0 <- function(a) 
  function(x) a[[1]]+x # + a[[2]]*x # + a[[3]]*x^2
warp1 <- function(a) 
  function(x) a[[1]]+ + a[[2]]*x # + a[[3]]*x^2

optim_fun <- optim_generator(warp0)

opt_res <- optim( 0, optim_fun, gr=NULL, clB, hessian = F, method='Brent', lower=-1000, upper=1000)
#Rprof()
#summaryRprof('work/profile')

startV <- c(0)
warp_par <- function(i){
  cspec <- readInd.TofH5(myTof, i)
  spspec <- semisparse_spec(cspec, lower=0, minlen=10, max_gap=30)
  sp2 <- sparse_remove_sat(spspec, sats)
  opt_res <- optim(startV, optim_fun, NULL, sp2, method='Brent', lower = -1000, upper=1000)
  #opt_res <- optim(startV, optim_fun, NULL, sp2, method='Nelder-Mead')
  list(index=i, opt=opt_res)
}

i_sel <- floor(seq(from =1 , to=myTof$indexhelp$N, length.out=20))
system.time( ww <- lapply(i_sel, warp_par))

(has_converged <- sapply(ww, function(x) x$opt$convergence))
# saturated peak make > 50% correlation
sapply(ww, function(x) x$opt$value)
pars <- sapply(ww, function(x) x$opt$par)
plot(pars) # kein shift!!!

plot(pars[1,])
plot(pars[2,])

# seem to have no drift
plot(pars[1,], pars[2,])

# pars 2 vs 3 klar linear -> 1 freiheitsgrad zu viel
plot(pars[3,])
plot(pars[3,], pars[2,])


