
tof.h5 <- 'testdata/Ac just breath after C (2014-10-23T11h34m53_#).h5'
fid <-H5Fopen(tof.h5)
tofblock <- raw_tofblock(fid)
indexhelp <- tof.indexhelp(tofblock)
spec0 <- read_spec_at(tofblock, indexhelp, 1)
spec1 <- read_spec_at(tofblock, indexhelp, 700)
spec2 <- read_spec_at(tofblock, indexhelp, 800)
spec3 <- read_spec_at(tofblock, indexhelp, 900)
spec4 <- read_spec_at(tofblock, indexhelp, 901)
plot(spec1-spec2,type='l')
