library(ptw)
library(dplyr)
library(parallel)
library(nloptr)

tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
myTof <- tof_h5(tof.h5)
aSpec <- read_spec_ind.tof_h5(myTof,10)
bSpec <- read_spec_ind.tof_h5(myTof,10000)
totalSpec <- sum_spec.tof_h5(myTof)

a2 <- aSpec/max(aSpec)
b2 <- bSpec/max(bSpec)
t2 <- totalSpec/max(totalSpec)

plot(t2, type='l')
lines(a2, col='red')

system.time(
ap <- ptw(t(t2),t(b2), optim.crit='WCC', trwdth=5, verbose=T)
)

#tof.h5 <- 'testdata/Ac just breath after C (2014-10-23T11h34m53_#).h5'
tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
fid <-H5Fopen(tof.h5)
tofblock <- raw_tofblock(fid)

# takes 30 s and 9 GB in process - maybe can be accelerated
# for the breath sample even 240 s = 4 min
system.time(
tof.spectra <- get.full.specblock(tofblock) 
)
#ass <- Matrix(tof.spectra, sparse=T)
#
#

# next idea make it iterative
totSigs <- colSums(tof.spectra)
tosum <- ( rowSums(tof.spectra))


# split matrix
rescale.prep <- list(inds = seq_along(tosum), sqind = sqrt(seq_along(tosum)))
re.index <- with(rescale.prep, function(ll)  ll[[1]] + ll[[2]]*sqind + ll[[3]]*inds)
ap.spec <- function(spc){ 
  cusp <- cumsum(spc)
  approxfun(x=seq_along(cusp), y=cusp, yleft=0, yright=cusp[length(cusp)])
}

prep.ref <- c(as.array(diff(tosum)),0)

ssRef = sqrt(sum(tosum^2))
nspec <- ncol(tof.spectra)
nsplits <- 7

dsplit <- ceiling(nspec/nsplits)
seqs <- seq(dsplit, nspec-1, by=dsplit)
ranges <- mapply(function(a,b) list(a:b),c(1,seqs+1), c(seqs,nspec))
specs <- lapply(ranges, function(x) rowSums(tof.spectra[,x]))
#data.frame(starts=c(1,seqs+1), ends=c(seqs,nspec))
pr <- 77700:78000
plot(pr,tosum[pr], type='l')
lines(pr,specs[[3]][pr], col='red')

aspec <- specs[[4]]

op_spec <- function(aspec){
  spec.fun <- ap.spec( aspec)
  warper <- function(ll){ diff(c(0,spec.fun(re.index(ll)))) }
  optfun <- function(ll){ 
        wspc <- warper(ll) 
        ssw <- sqrt(sum(wspc^2))
        -crossprod(tosum, wspc)/ssw/ssRef 
  }

  # see: http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
  #op.ra <- nloptr(c(0,0,1), optfun, 
                                  #opts = list("algorithm"="NLOPT_LN_COBYLA", 
                                                                #ftol_rel=1e-6, xtol_rel=1e-6, maxeval=1000),
                                  #lb=c(-1e6,-1e4,0), ub=c(1e6,1e4,2))

  op.rao <- optim(c(0,0,1), optfun) 
  ws <- warper(op.rao$par)

#pr <- 230000:240000
#plot(pr,tosum[pr], type='l')
#lines(pr,ws[pr], col='blue')

  list(parameter = op.rao$par,
           warped = ws)
}

system.time(
w1 <- lapply(specs, op_spec)
)

tmp <- do.call(rbind,  lapply(w1, function(x) with(x,warped)) )

s2 <- apply(tmp,2,sum)
ss2 <- smooth(s2)
#pr <- 178000:178500
pr <- 230000:240000
pr <- 222000:222500
pr <- 236400:236700
plot(pr,s2[pr], type='l')
lines(pr,ss2[pr], col='blue')
lines(pr,tosum[pr], col='red')
lines(pr,tmp[7,pr], col='cyan')
#lines(pr,tmp[4,pr], col='cyan')

wf <- do.call(rbind,  lapply(w1, function(x) with(x,parameter)) )
plot(wf[,1])
plot(wf[,2])
plot(wf[,3])


#global
#op.ra <- nloptr(c(0,1,0), optfun, 
                                #opts = list("algorithm"="NLOPT_GN_DIRECT_L", ftol_rel=1e-4),
                                #lb=c(-1e6,0,0), ub=c(1e6,2,1e-3))

# ---NLOPT_GN_DIRECT_L
#split.at <- floor(ncol(tof.spectra)/2)
#left.spec <- tof.spectra[,1:split.at]
#right.spec <- tof.spectra[,-(1:split.at)]
#lefts <- rowSums(left.spec)
#rights <- rowSums(right.spec)

#pr <- 178000:178500
#plot(pr,tosum[pr], type='l')
#lines(pr,lefts[pr], col='red')
#lines(pr,rights[pr], col='blue')

left.fun <- ap.spec(lefts)
right.fun <- ap.spec(rights)

optfun <- function(ifun,ll){
  i2 <- re.index(ll)
  s2 <- ifun(i2)
  crossprod(prep.ref,s2)
}

o1 <- optim(c(0,1,0), function(x) optfun(left.fun,x))
o2 <- optim(c(0,1,0), function(x) optfun(right.fun,x))

left2 <- diff(c(0, left.fun(re.index(o1$par))))
right2 <- diff(c(0, right.fun(re.index(o2$par))))

pr <- 178000:178500
plot(pr,tosum[pr], type='l')
lines(pr,left2[pr], col='red')
lines(pr,lefts[pr], col='blue')

plot(pr,tosum[pr], type='l')
lines(pr, right2[pr], col='red')
lines(pr, rights[pr], col='blue')


pr <- 171000:171500
plot(pr,tosum[pr], type='l')
lines(pr, (left2+right2)[pr], col='red')

# old
#

n<- 1000
sumspec <- read.sum.spec(fid)
cumspec <- cumsum(sumspec)

master.app <- approxfun(x=seq_along(cumspec), y=cumspec, yleft=0, yright=cumspec[length(cumspec)])


ss <- tof.spectra[,234,drop=F]
reference <- Matrix(sumspec, sparse=T)

system.time(
test <- ptw(t(sumspec), t(tof.spectra[,1:10]), optim.crit='RMS', smooth.param=0)
)
plot(test)

pr = seq(200000, 220000)
ms <- max(reference[pr])
mc <- max(test$warped.sample[,pr])
ws <- test$warped.sample[,pr]
wo <- test$sample[,pr]
plot(pr,reference[pr]/ms, type='l', lwd=4)
i <- 5
lines(pr,ws[i,], col='red')
lines(pr,wo[i,], col='blue')

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
