library(Rcpp)
library(dplyr)
library(splines)
library(ggplot2)
library(OrgMassSpecR)
library(MASS)
sourceCpp(file="src/peakFind.cpp")
electron.mass <- 1/1820

r.path <- file.path('R')
sl <- file.path(r.path,list.files(r.path, pattern='.*\\.R$'))
lapply(sl, source)

# consider: 
#Improved Calibration of Time-of-Flight Mass Spectra by Simplex Optimization of Electrostatic Ion Calculations

#Noah P. Christian ,† Randy J. Arnold ,‡ and James P. Reilly *
#Department of Chemistry, Indiana University, Bloomington, Indiana 47405
#Anal. Chem., 2000, 72 (14), pp 3327–3337
#DOI: 10.1021/ac991500h

tof.h5 <- file.path('testdata/',"2015.07.17-10h40m34 Ethanol deurated Karl .h5")

h5.holder <- function(tof.h5){
  fid <-H5Fopen(tof.h5)
  tofblock <- get.raw.tofblock(fid)
  indexhelp <- tof.indexhelp(tofblock)
  cr <- make.curr.tofreader(tofblock, indexhelp)

  list(close = function() H5close(),
	   readIndex = cr, 
	   sampleScans = function(n=100) {
		 with(indexhelp, 
		 	  if(N > 2*n) { floor(seq(from=1,to=N, length.out=n))}
			  else {	seq(from=1,to=N)}
					)},
	   read.scan = function (i) read.spec.ind(tofblock, indexhelp, i)
	   )
	  
}

curr.h5 <- h5.holder(tof.h5)
#mc.table <- mass.calib.tof(tofblock, indexhelp)

kn <- curr.h5$sampleScans(50)

locate.ions.block <- function(df.ion,curr.h5,
  preliminary.coeff=list(intercept = 900,square_mass = 17650 ),
  n.samp=100 ) {
  scl <- curr.h5$sampleScans(n.samp)
  pos.block <- expand.grid(scan=scl , ion=df.ion$ions)
  pp <- apply(pos.block, 1,function (dd) { 
				curr.spec.line <- curr.h5$read.scan( dd[['scan']])
				pos <- find_ion_tof(dd[['ion']], preliminary.coeff, curr.spec.line)
  })

  pos.block$pos <- pp
  pos.block <- pos.block[!is.na(pos.block$pos),]
  
  pos.block
}

# TODO: Massen suchen die *eindeutig* sind für die Kalibration

pos.block <- locate.ions.block(df.ion,curr.h5, 
							   preliminary.coeff=list(intercept=1850, square_mass=17500))



pos.add.data <- function(pos.block) 
  within(pos.block,{
		 sq.mass <- sqrt(ion)
		 sq3.mass <- sq.mass*ion
		 sq5.mass <- sq3.mass*ion
  })

pos.block <- pos.add.data(pos.block)

nn <- 3
kn <- floor(seq(1,to=indexhelp$N, length.out=nn+2))[2:nn+1]
#mod <- rm(pos ~ (ion + sq.mass + sq3.mass + sq5.mass) * ns(scan,knots=kn), pos.block )
rmod <- rlm(pos ~ ion + sq3.mass + sq5.mass + sq.mass * ns(scan,knots=kn), pos.block, method='MM')
#pl <- predict(mod, newdata=data.frame(sq.mass=sqrt(47.05), scan=1:indexhelp$N))
#plot(pl, type='l')


pl <- predict(rmod, newdata=pos.block)
pos.block$pred <- pl
ggplot(pos.block) + geom_line(mapping=aes(x=scan,y=pos,color=as.factor(ion))) + geom_line(mapping=aes(x=scan,y=pred,color=as.factor(ion))) + facet_wrap(~ion, scales='free_y')


ggplot(pos.block) + geom_point(mapping=aes(x=pos,y=pred,color=as.factor(ion))) + facet_wrap(~ion, scales='free')

  #data.frame(scan = (1:indexhelp$N) ) %>% rowwise() %>% do( {
    #curr.spec.line <- read.spec.ind(tofblock, indexhelp, .$scan)
    #mc <- mass.calib.coeff.single(ions, preliminary.coeff=pars.approx, curr.spec.line)
    #data.frame(scan=.$scan, intercept=mc[['intercept']], square_mass=mc[['square_mass']])


ind <- 50
masslist <- seq(46.7,47.3, by=0.001)
pl <- predict(rmod, newdata=read.block)
sc <- curr.h5$read.scan(ind)
rs <- readScales(sc,pl)

pd <- data.frame(mass=masslist, cr=rs)

ggplot(pd) + geom_line(mapping=aes(x=mass, y=cr)) + scale_y_log10()


# plot ethanol
masslist <- seq(46.9,47.1, by=0.001)

x <- filter(read.block, scan==10)



read.scan.data.sync <- function(curr.h5, rmod, masslist,n.scans=50){
  with(curr.h5,{
		 streaming.scan <- function(x){
		   sc <- read.scan( x[[1,'scan']])
		   rs <- readScales(sc,x$inds)
		   x$cr <- rs
		   return(x)
		 }
  		 scan.ind <- sampleScans(n.scans)
		 dfm <- expand.grid(ion=masslist, scan=scan.ind)
		 read.block <- pos.add.data(dfm)
		 read.block$inds <- predict(rmod, newdata=read.block)
		 read.scans <- read.block %>% group_by(scan) %>% do( streaming.scan(.))
  })
}
dark.side.plot <- function(read.scans) {
  nsc <- max(read.scans$scan)
  vs <- max(read.scans$cr, na.rm=TRUE)/3/nsc
  hs <- (max(read.scans$ion)- min(read.scans$ion))/nsc/5
  pd <- read.scans %>% mutate(x.plot=ion + scan*hs, y.plot=cr + scan*vs)
  ggplot(pd) + 
	geom_line(mapping=aes(x=x.plot, y=y.plot, group=as.factor(scan), 
  							color=scan), alpha=0.3) + 
  	  theme_bw() + 
  	  	scale_y_continuous('counts') + scale_x_continuous('mass')+
  	  	  scale_color_gradient(guide=FALSE, low='red', high='blue')
}

#ggplot(pd) + geom_line(mapping=aes(x=x.plot, y=y.plot, color=as.factor(scan)), alpha=0.3) + theme_bw()
masslist <- seq(46.95,47.05, length.out=200)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist)
dark.side.plot(rsd)

# etoh D6
masslist <- seq(52.95,53.1, length.out=200)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist)
dark.side.plot(rsd)
MonoisotopicMass(formula=list(C=5,H=8),charge=1)
masslist <- seq(46,54, by=0.001)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist, n.scans=20)
dark.side.plot(rsd)
masslist <- seq(48.95,49.05, length.out=200)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist)
dark.side.plot(rsd)

# acetald D4
masslist <- seq(49.95,50.05, length.out=200)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist)
dark.side.plot(rsd)

# acetat normal
MonoisotopicMass(formula=list(C=2,H=4,O=2),charge=1)
masslist <- seq(61.0,61.1, length.out=200)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist,n.scans=40)
dark.side.plot(rsd)
# acetat D4
masslist <- seq(65.0,65.1, length.out=200)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist)
dark.side.plot(rsd)
# acetat D3
masslist <- seq(63.9,64.2, by=0.001)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist, n.scans=20)
dark.side.plot(rsd)

# isopren 
MonoisotopicMass(formula=list(C=5,H=8),charge=1)
masslist <- seq(68.9,69.2, by=0.001)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist, n.scans=20)
dark.side.plot(rsd)
MonoisotopicMass(formula=list(C=5,H=8),charge=1)
masslist <- seq(68,75, by=0.001)
rsd <- read.scan.data.sync(curr.h5, rmod, masslist, n.scans=20)
dark.side.plot(rsd)
