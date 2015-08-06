H5close()
tof.h5 <- file.path('testdata/Ac just breath after D (2014-10-23T13h54m32_#).h5')

n <- new("TofH5", file.name=tof.h5)
n2 <- open(n)
n3 <- close(n2)


n2 <- open(n)
data(sample_ions)


read.scan(n2,1)
