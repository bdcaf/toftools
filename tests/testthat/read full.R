library(dplyr)
library(tidyr)

tof.h5 <- file.path('data',"uncal.h5")
h5ls(tof.h5)

fm <- h5read(tof.h5,'FullSpectra/TofData')
class(fm)

flat.tof <- read.tof.fullspectra(file.path('data',"uncal.h5"))
mc.table <- mass.calib.tof(flat.tof)
smooth.mc <- smooth.mass.cal(mc.table)
