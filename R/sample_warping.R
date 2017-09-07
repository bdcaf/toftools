source('R/tofReader.R')
source('R/tofSparse.R')
library(zoo)
#tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
tof.h5 <- 'testdata/2017.06.22-11h13m23s ca valve open.h5'
myTof <- tof_h5(tof.h5)
aSpec <- read_spec_ind.tof_h5(myTof,10)
full.wave <- aSpec
which(diff(full.wave) < -1000)
plot(full.wave, type='l')
plot(diff(full.wave), type='l')
plot((full.wave), type='l', xlim=c(79000,81000))
plot(diff(full.wave), type='l', xlim=c(9.96e4, 9.98e4))
system.time({
spsp <- semisparse_spec(full.wave, lower=0, minlen=10, max_gap=30)
})
#system.time({
#simplified <- simplify_sparse(spsp, max_gap=50L)
#})

warp0 <- function(a) 
  function(x) a[[1]]+x # + a[[2]]*x # + a[[3]]*x^2

totalSpec <- sum_spec.tof_h5(myTof)
full.wave2 <- totalSpec / myTof$indexhelp$N
plot(full.wave2, type='l')
which(diff(full.wave2) < -300)
plot(diff(full.wave2), type='l')
plot((full.wave), type='l', xlim=80007+c(-100,+200))
plot((full.wave), type='l', xlim=99707+c(-100,+200))
plot((full.wave), type='l', xlim=102855+c(-100,+200))

refSpec <- smooth(totalSpec, twiceit=T)
sats <- find_saturated(refSpec, N=myTof$indexhelp$N)
ref2 <- dense_remove_sat(refSpec, sats)

#refSpec <- rollmean(totalSpec, 10)
#refSpec <- totalSpec
#refSpec <- totalSpec
refEn <- sum(ref2^2)

normSpec <- ref2/sqrt(refEn)

assure_range <- function(f)
   function(x) pmax.int(pmin.int(f(x),length(normSpec)),1)



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
cspec <- read_spec_ind.tof_h5(myTof, i)
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

startV <- c(0)
warp_par <- function(i){
  cspec <- read_spec_ind.tof_h5(myTof, i)
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


