
library(rhdf5)
file.path('data',"Ac just exhale (2014-10-23T10h43m42_#).h5")

fid <-H5Fopen(tof.h5)
fid
at <- H5Aopen(fid,"NbrWaveforms")
waveforms <- H5Aread(at)


gr <- H5Gopen(fid,"FullSpectra")
tof.data <- H5Dopen(gr,"TofData")
tof.data

h5spaceFile <- H5Dget_space(tof.data)
dims <- H5Sget_simple_extent_dims(h5spaceFile)
index = NULL

calc.indices <- expand.grid(buf=1:dims$size[[3]], write=1:dims$size[[4]])


pos <- calc.indices[235,]

# h5spaceMem <- H5Screate_simple(size)

h5spaceMem <- H5Screate_simple(dims$size[[1]])
# H5Sselect_hyperslab(h5spaceFile, 
#                     start = c(1,1,pos$buf,pos$write), stride = c(1,1,1,1), count = c(dims$size[[1]],1,1,1), 
#                     block = NULL)
size <-  H5Sselect_hyperslab(h5spaceFile, start = c(1,1,pos$buf,pos$write), count = c(dims$size[[1]],1,1,1) )                          
single.tof <- H5Dread(h5dataset = tof.data, h5spaceFile = h5spaceFile, 
               h5spaceMem = h5spaceMem )
plot(single.tof, type='l')


selected <- H5Sselect_index(tds,index=list(1:4,1,1,1))
selected
#398999
target <- H5Screate_simple(dims <- list(398999), maxdims  <- list(398999))


r <- H5Dread(tof.data,h5spaceFile = tds, h5spaceMem = target)
?H5Dread
