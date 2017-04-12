source('R/spec_tools.R')
source('R/tofReader.R')
source('R/massCalTools.R')
library(dtw)
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

inds <- seq_len(indexhelp$N) 
aspec <- function(n) read.spec.ind(tofblock, indexhelp, n)
rf <- function(x,n) x + aspec(n)
#ss <- Reduce(rf, inds, 0, accumulate=F)
ss <- read.sum.spec(fid)
plot(ss, type='l')

m2 <- masscal_legacy_smooth(fid, indexhelp)


make_massaxis <- function( resolution=1e4,
						   mass.low=10,
						   mass.high=200){
  log_mass_axis <- seq(from=log10(mass.low), 
					   to = log10(mass.high),
					   length.out=((mass.high-mass.low)*resolution)
					   )
  mass_axis <- 10^log_mass_axis
}

calib <- m2[1,]
spec <- aspec(1)
ma <- with(calib, (seq_along(spec)-b)^2/a^2)
plot(ma, spec, type='l')

spec2mav <- function(spec, calib, mass_axis){
  sqrt(mass_axis)
}

masscal_c1 <- function(fid, indexhelp, 
					   resolution=1e4,
					   mass.low=10,
					   mass.high=200){
  mass_axis <- make_mass_axis(resolution=resolution, mass.low=mass.low, mass.high=mass.high)
}


# mass calibration jumps
legacy_masscal <- read.mass.cals(fid)
plot(legacy_masscal[1,], type='l')
plot(legacy_masscal[2,], type='l')


# mass calibration a smooth one makes maybe sense
df <- as_data_frame(t(legacy_masscal)) %>% mutate(write = 1:nrow(.))
a_smooth <- with(df, smooth.spline(write, V1))
b_smooth <- with(df, smooth.spline(write, V2))


r <- with(indexhelp, 1/dims[[3]])
ci <- indexhelp$calc.indices %>% 
  mutate(x = write - 0.5 + r*buf,
  		 a = predict(a_smooth, x=x)$y,
  		 b = predict(b_smooth, x=x)$y
  		 )

par(mfrow=c(2,1))

plot(df$write, df$V1, type='p')
lines(ci$x, ci$a)

plot(df$write, df$V2, type='p')
lines(ci$x, ci$b)


left_join(indexhelp$calc.indices, df, by='write')
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

s1 <- sparse_spec(spec0)
s2 <- sparse_spec(spec4)
	
## algorithm 2

len1 <- sum(s1$rlen)
len2 <- sum(s2$rlen)

df <- function(v,w) abs(v-w)
phif <- function(s,t) 


mm <- matrix(data = Inf, nrow=len1, ncol=len2)
mm[1,1] <- 0


v1 <- s1$ser[[2]]
v2 <- s2$ser[[4]]

tw1 <- dtwDist(v1,v2)

# idee: benutze convolve für lokale peak verschiebungen
cc <- convolve(v1,v2, type='open')
c0 <- convolve(v1,v1, type='open')

to_shift <- which.max(cc) - length(v1)


peak1v <- s1$ser[[2]]
peak2v <- s2$ser[[4]]
costP12 <- function(peak1v,peak2v){
  cc <- convolve(peak1v,peak2v, type='open')
  mm <- max(cc)
  cost <- 1 - cc/mm
}


rec_dist <- function(distm){
  if (is.null(distm) || is.null(dim(distm))) return(NULL)

  ind <- which.min(distm)
  if (length(ind) < 1) return(NULL)
  subs <- arrayInd(ind[[1]], dim(distm))

  bef_mat <- distm[1:subs[[1]]-1,1:subs[[2]]-1]
  aft_mat <- distm[(subs[[1]]+1): nrow(distm),(subs[[2]]+1): ncol(distm)]

  #nds = which(distm == min(distm), arr.ind=TRUE)
  #o0 <- data.frame(nds)
  #distm[o0$row,] <- Inf
  #distm[,o0$col] <- Inf

  #o1 <- rec_dist(distm)

  bef_res <- rec_dist(bef_mat)
  aft_res <- rec_dist(aft_mat)
  if (!is.null(aft_res)){
	aft_res[,1] <- aft_res[,1]+subs[[1]]
	aft_res[,2] <- aft_res[,2]+subs[[2]]
  }
  rbind(bef_res, subs) %>% rbind(aft_res)
}


plot(raw_peak_comb)

pair_peaks <- function(s1,s2, maxDist = 1e3){
  dm <- expand.grid(s1 = s1$glob_max, s2=s2$glob_max) %>%
	mutate(dist = (s1 - s2)^2, d2 = ifelse(dist>maxDist, Inf, dist))
  distm <- matrix(data = dm$d2, nrow = nrow(s1), ncol = nrow(s2))

  raw_peak_comb <- rec_dist(distm)

  colnames(s1) <- paste(colnames(s1),'1')
  colnames(s2) <- paste(colnames(s2),'2')
  paired <- with(raw_peak_comb, as_data_frame(
						  cbind(s1[row,], 
								s2[col,])))
}


raw_pm <- pair_peaks(s1, s2)
