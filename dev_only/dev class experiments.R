source('R/tofReader.R')
source('R/TofMeasurementClass.R')
source('R/TofGetScan.R')
source('R/TofGetJustScan.R')
source('R/massCalFunctions.R')

tof.h5 <- file.path('testdata',"Ac just exhale (2014-10-23T10h43m42_#).h5")


tm <- createTofMeasurement(file=tof.h5)
sc <- getJustScan(tm,1000)

plot(sc, type='l')

spectra <- transformIntensity(sc, method="sqrt")

spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=3)

plot(sc, xlim=c(68,70), ylim=c(0,10))
lines(spectra, col='red')

# baseline <- estimateBaseline(spectra, method="SNIP", iterations=100)
# 
# plot(baseline)
#
# is all 0!


# library(xcms)
# detach(name="package:xcms", unload = TRUE)

peaks <- detectPeaks(sc, method="SuperSmoother", SNR=60, halfWindowSize=6)
plot(sc, xlim=c(59,59.1))
points(peaks)

cbind(peaks@mass, peaks@intensity, peaks@snr)
# line(baseline)

range <- 1:tm@.indexHelp$N

get.maxima <- function(i,tm, win=15) {
  sc <- getJustScan(tm,i)
  logical.max <- MALDIquant:::.localMaxima(sc,win) 
  # use halfWindowSize of 6
  pos <- which(logical.max)
  good <- sc[pos] > 10
  # use 10 counts to identify *landmarks*
  list(pos[good])
}

get.maxima(2000,tm)

range=1:100
all.max <- lapply(range, function(x) get.maxima(x,tm))

# need to bin peaks

sp <- sort(unlist(all.max), method="quick")
dsp <- diff(sp)
candidates <- dsp>15
cut.cand <- c(0,sp[candidates]+1,Inf)
which(candidates)
sp[1:20]
grouped <- cut(sp,cut.cand)
library(dplyr)
library(tidyr)
an <- data.frame(grouped,sp)

a <- tm@metaDataScan[['MassCalibration a']]
b <- tm@metaDataScan[['MassCalibration b']]
# mv <- .idx2mass(1:length(spec),a,b)
mass.suggestion <- an %>% group_by(grouped) %>% 
  summarise(mean=mean(sp), median=median(sp), n=n()) %>% filter(n>70) %>%
  mutate(mass.pre = .idx2mass(mean,a,b))