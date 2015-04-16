source('dev_only/dev class.R')
source('dev_only/TofGetScan.R')
source('R/tofReader.R')

tof.h5 <- file.path('testdata',"Ac just exhale (2014-10-23T10h43m42_#).h5")


tm <- createTofMeasurement(file=tof.h5)
sc <- getScan(tm,100)

plot(sc, type='l')
