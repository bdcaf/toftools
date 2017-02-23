library(dtw)

tof.h5 <- 'testdata/Ac just breath after C (2014-10-23T11h34m53_#).h5'
fid <-H5Fopen(tof.h5)
tofblock <- get.raw.tofblock(fid)
indexhelp <- tof.indexhelp(tofblock)
spec0 <- read.spec.ind(tofblock, indexhelp, 1)
spec1 <- read.spec.ind(tofblock, indexhelp, 700)
spec2 <- read.spec.ind(tofblock, indexhelp, 800)
spec3 <- read.spec.ind(tofblock, indexhelp, 900)
spec4 <- read.spec.ind(tofblock, indexhelp, 901)
plot(spec1-spec2,type='l')

range = 136500:136600

spectmp <- lapply(1:10, function(x) read.spec.ind(tofblock, indexhelp, x))
specmat <- do.call(rbind, spectmp)

candval <- apply(specmat,1,max)


plot(seq_along(spec1),spec1,type='l', xlim=c(1.365e5,1.366e5), ylim=c(0,40))
lines(seq_along(spec2),spec2, col='red')

w1 <- spec1[range]
w2 <- spec2[range]
w3 <- spec3[range]
w4 <- spec4[range]
alignment <- dtw(w3, w4, 
				 dist.method = 'correlation',
				 window.type = "sakoechiba", 
				 window.size=5, 
				 keep=T)
plot(alignment, type="two")
plot(alignment)

plot(spec0-spec1,type='l')
#plot(spec1,spec2,type='p')

## from example
## A noisy sine wave as query
idx<-seq(0, 100, len=500);

## A cosine is for reference; sin and cos are offset by 25 samples
reference <- 2*dnorm(idx, mean=10, sd=1) +
			 dnorm(idx, mean=30, sd=2) +
			 0.5*dnorm(idx, mean=70, sd=3)

query     <- 0.9*dnorm(idx, mean=12, sd=1.2) +
			 3*dnorm(idx, mean=34, sd=2.25) +
			 2*dnorm(idx, mean=68, sd=2.2) +
			 0.6*dnorm(idx, mean=24, sd=2.2) +
			 0.1*rpois(length(idx), 0.2)
		  
plot(reference); lines(query,col="blue"); 

## Find the best match
alignment<-dtw(query,reference, keep=T,
			   open.begin=T,
			   open.end=T,
			   step.pattern=asymmetric)

## Display the mapping, AKA warping function - may be multiple-valued
## Equivalent to: plot(alignment,type="alignment")
plot(alignment$index1,alignment$index2,main="Warping function");
plot(alignment, type='two')

## Confirm: 25 samples off-diagonal alignment
lines(1:100-25,col="red")

summary(pr_DB)
#alignment<-dtw(query,reference, dist.method='cosine', keep=T);
#alignment<-dtw(query,reference, dist.method='Braun-Banquet', keep=T);
alignment<-dtw(query,reference, 
				keep=T,
				dist.method='cosine',
				open.begin=F,
				open.end=F,
				step.pattern= symmetric2,
				window.type= 'sakoechiba',
				window.size = 50)
plot(alignment, type='two')
plot(alignment)

# blind DTW scheint schlechte Lösung zu sein
# vielleicht mit Massenkalibration so dass overlap maximal ist.
