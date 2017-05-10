source('R/spec_tools.R')
source('R/tofReader.R')
source('R/massCalTools.R')
library(ptw)
library(Matrix)
library(dplyr)
library(parallel)

tof.h5 <- 'testdata/Ac just breath after C (2014-10-23T11h34m53_#).h5'
fid <-H5Fopen(tof.h5)
tofblock <- get.raw.tofblock(fid)

# takes 30s and 9GB in process - maybe can be accelerated
system.time(
tof.spectra <- get.full.specblock(tofblock) 
)
#ass <- Matrix(tof.spectra, sparse=T)

n<- 1000
sumspec <- read.sum.spec(fid)
cumspec <- cumsum(sumspec)

master.app <- approxfun(x=seq_along(cumspec), y=cumspec, yleft=0, yright=cumspec[length(cumspec)])

ifun <- function(ll, inds)  ll[[1]] + ll[[2]]*inds + ll[[3]]*inds^2

ss <- tof.spectra[,234,drop=F]
reference <- Matrix(sumspec, sparse=T)

system.time(
test <- ptw(t(sumspec), t(tof.spectra[,1:10]))
)
plot(test)

pr = seq(152000, 160000)
ms <- max(sumspec[pr])
mc <- max(ss[pr])
plot(pr,sumspec[pr]/ms, type='l', lwd=4)
lines(pr,ss[pr]/mc, col='red')

targets <- 0:length(cumspec)
newind <- function(ll) ifun(ll, targets)
resampe <- function(ll){
  i2 <- newind(ll)
  diff(master.app(i2))
}
opt_shift <- function(n){ 
  sa <- as(tof.spectra[,n,drop=F], 'sparseVector')
  optim.fun <- function(ll){
	-crossprod(sa*resampe(ll))[1]
  }

  best <- optim(c(0.,1,0.), optim.fun, control=list())
  with(best, par)
}

check_inds <- seq(from=1, to=indexhelp$N, length.out=30)
system.time(
  v <- mclapply(check_inds, opt_shift, mc.cores=8L, mc.preschedule=F)
)

vs <- sapply(v, as.numeric)
#vs <- do.call(rbind, v2)
par(mfrow=c(3,1))
plot(vs[1,] ,type='l')
plot(vs[2,] ,type='l')
plot(vs[3,] ,type='l')
