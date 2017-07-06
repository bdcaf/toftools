source('R/tofReader.R')
tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
myTof <- tofH5(tof.h5)
aSpec <- readInd.TofH5(myTof,10)
full.wave <- aSpec
system.time({
spsp <- sparse_spec(full.wave, lower=0, minlen=10, max_gap=30)
})
#system.time({
#simplified <- simplify_sparse(spsp, max_gap=50L)
#})

spt <- sparseTof(aSpec)
totalSpec <- sumSpec.TofH5(myTof) # sparsity makes no sense here
refE <- sum(totalSpec^2)
system.time(
cor.semisparse.full(spsp, totalSpec, refE)
)

warp0 <- function(a) 
  function(x) a[[1]] + a[[2]]*x

a <- c(10,0.99)
refSpec <- totalSpec
refEn <- sum(refSpec^2)
optim_fun <- function(a){
  wf <- warp0(a)
  warped <- warp_spec(spsp, wf)
  -cor.semisparse.full(warped, totalSpec, refEn)
}


r <- optim_fun(a)
system.time(
opt_res <- optim(c(0,1), optim_fun)
)


opt_res

