source('R/tofReader.R')
source('R/TofMeasurementClass.R')
source('R/TofGetScan.R')
source('R/massCalFunctions.R')

tof.h5 <- file.path('testdata',"Ac just exhale (2014-10-23T10h43m42_#).h5")


tm <- createTofMeasurement(file=tof.h5)
sc <- getScan(tm,100)

plot(sc)

# spectra <- transformIntensity(sc, method="sqrt")

spectra <- smoothIntensity(sc, method="SavitzkyGolay", halfWindowSize=3)

plot(spectra)


# baseline <- estimateBaseline(spectra, method="SNIP", iterations=100)
# 
# plot(baseline)
#
# is all 0!
