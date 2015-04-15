
tof.h5 <- file.path('data',"Ac just exhale (2014-10-23T10h43m42_#).h5")

ex <- make.tof.h5(tof.name = tof.h5)
close(ex)          

ex
ex <- open(ex)
read.waveforms(ex)
ex <- close(ex)

ex          
