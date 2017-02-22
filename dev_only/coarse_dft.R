source('R/tofReader.R')
library(dplyr)
# load the data
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

# trying to follow
# 1. Dupont, M. & Marteau, P.-F. in Lecture Notes in Computer Science
# (including subseries Lecture Notes in Artificial Intelligence and
# Lecture Notes in Bioinformatics) (eds. Douzal-Chouakria, A., Vilar,
# J. A. & Marteau, P.-F.) 9785, 157–172 (Springer International
# Publishing, 2016).
#
# sparse time series:
# 
full.wave <- as.vector(spec1)

sparse_spec <- function(full.wave){
  idx <- full.wave > 0
  renc <- rle((idx))
  startp <- cumsum( c(0, as.numeric(renc$lengths)))

  extract_ts <- function(starter, len)
	list(start = starter, 
		 rlen = len,
		 signal = full.wave[starter + seq_len(len)])

  dv <- data.frame(start=startp[seq_along(renc$lengths)], 
				   rlen = renc$lengths, 
				   greater0 = renc$values) %>%
				filter(rlen >3, greater0) 

  mapply(extract_ts, dv$start, dv$rlen)
}
	
## algorithm 2
