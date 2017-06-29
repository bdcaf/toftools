source('R/spec_tools.R')
source('R/tofReader.R')
source('R/massCalTools.R')
library(dplyr)
library(PeakSegDP)

#tof.h5 <- 'testdata/Ac just breath after C (2014-10-23T11h34m53_#).h5'
tof.h5 <- 'testdata/2017.02.15-15h22m12s D6-EtOHbreathclemens.h5'
myTof <- tofH5(tof.h5)
aSpec <- readInd.TofH5(myTof,10)
totalSpec <- sumSpec.TofH5(myTof)

plot(totalSpec/max(totalSpec), type='l')
lines(aSpec/max(aSpec), col='red')


system.time(
tmp <- cDPA(as.vector(aSpec), maxSegments = 20)
)

# takes 30 s and 9 GB in process - maybe can be accelerated
# for the breath sample even 240 s = 4 min
system.time(
tof.spectra <- get.full.specblock(tofblock) 
)


