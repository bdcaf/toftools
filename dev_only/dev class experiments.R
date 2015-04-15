tof.h5 <- file.path('testdata',"Ac just exhale (2014-10-23T10h43m42_#).h5")


tm <- createTofMeasurement(file=tof.h5)
sc <- getScan(tm)

plot(sc, type='l')
