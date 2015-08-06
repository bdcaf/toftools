H5close()
tof.h5 <- file.path('testdata/Ac just exhale (2014-10-23T10h43m42_#).h5')

n <- new("TofH5", file.name=tof.h5)
n2 <- open(n)
n3 <- close(n2)


n2 <- open(n)
data(sample_ions)
ions <- sample_ions


sc.test<- read.raw.scan(n2,100)


####
imc <- tofTools:::init.mass.calib(n2)
x <- with(imc, ((seq_along(sc.test)-intercept)/square_mass)^2)

ig <- x > 30 & x <40

plot(x[ig], sc.test[ig], type='l')


mass.calibration <- mass.calib(n2, sample_ions)

rb <- read.bins(n2, mass.calibration, index=1:10, masses=20:50)
