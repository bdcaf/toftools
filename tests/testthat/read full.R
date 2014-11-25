library(dplyr)
library(tidyr)

list.files('data')
tof.h5 <- file.path('data',"uncal.h5")
h5ls(tof.h5)

fm <- h5read(tof.h5,'FullSpectra/TofData')
class(fm)

tof.h5 <- file.path('data',"Ac just exhale (2014-10-23T10h43m42_#).h5")
h5ls(tof.h5)

flat.tof <- read.tof.fullspectra(file.path('data',"Ac just exhale (2014-10-23T10h43m42_#).h5"))
mc.table <- mass.calib.tof(flat.tof)
smooth.mc <- smooth.mass.cal(mc.table)

mass.list <- c(21.022, 59.04914, 69.06988)
mr1 <- lapply(mass.list, mass.to.range)

mass.range <- mr1[[2]]
c.ind <- 5
curr.mass.cal <- smooth.mc[c.ind,]

c.spec <- flat.tof[,c.ind]
c.range <- mass.range2tof.range(curr.mass.cal, mass.range)
ran <- c.range[['lower']]:c.range[['upper']]
plot(c.spec[ran], type='l')


one.peak.integrate(curr.mass.cal = c.mc, mass.range = cm, c.spec)

pa <- apply(smooth.mc,1, function(x){  
  c.spec <- flat.tof[,x['ind']]
  one.peak.integrate(curr.mass.cal = x, mass.range = cm, c.spec)
})

plot(pa, type='l')
