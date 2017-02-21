# [1] E. Moskovets and B. L. Karger, “Mass calibration of a matrix-assisted laser desorption/ionization time-of-flight mass spectrometer including the rise time of the delayed extraction pulse,” Rapid Commun. Mass Spectrom., vol. 17, no. 3, pp. 229–237, 2003.


enrich.ib <- function(ion.block){
  within(ion.block,{
    sq <- sqrt(ion)
    isq <- 1/sq
  })  
}

tof.h5 <- file.path('testdata/Ac just exhale (2014-10-23T10h43m42_#).h5')
data(sample_ions)
ions <- sample_ions

object <- new("TofH5", file.name=tof.h5)
object <- open(object)

scl <- sample.scans(object, n.samp=100)

ion.block <- find.ion.block(object,ions)


ion.block <- enrich.ib(ion.block)
head(ion.block)
num.knots <- 2
kn <- floor(seq(1,to=object@.indexhelp$N, length.out=num.knots+2))[2:num.knots+1]
fit.mos <- rlm(pos ~ (isq + ion + sq) * ns(scan,knots=kn), ion.block, method='MM')


extr.data <- expand.grid(scan= 1:10, ion=seq(47,47.1,length.out=100))
extr.data <- enrich.ib(extr.data)
extr.data$idx <- predict(fit.mos, newdata = extr.data)

read.prepared.bins(object,targets=extr.data)
