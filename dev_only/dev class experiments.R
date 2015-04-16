source('R/tofReader.R')
source('R/TofMeasurementClass.R')
source('R/TofGetScan.R')
source('R/massCalFunctions.R')

tof.h5 <- file.path('testdata',"Ac just exhale (2014-10-23T10h43m42_#).h5")


tm <- createTofMeasurement(file=tof.h5)
sc <- getScan(tm,100)

plot(sc)

spectra <- transformIntensity(sc, method="sqrt")

spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=30)

plot(spectra, xlim=c(58,60))


# baseline <- estimateBaseline(spectra, method="SNIP", iterations=100)
# 
# plot(baseline)
#
# is all 0!


# library(xcms)
# detach(name="package:xcms", unload = TRUE)

peaks <- detectPeaks(spectra, method="SuperSmoother", SNR=5, halfWindowSize=50)
plot(spectra, xlim=c(50,60))
points(peaks)

cbind(peaks@mass, peaks@intensity, peaks@snr)
# line(baseline)

