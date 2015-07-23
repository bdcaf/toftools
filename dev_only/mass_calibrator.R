library(splines)
library(ggplot2)
library(OrgMassSpecR)
library(MASS)
electron.mass <- 1/1820

# consider: 
#Improved Calibration of Time-of-Flight Mass Spectra by Simplex Optimization of Electrostatic Ion Calculations

#Noah P. Christian ,† Randy J. Arnold ,‡ and James P. Reilly *
#Department of Chemistry, Indiana University, Bloomington, Indiana 47405
#Anal. Chem., 2000, 72 (14), pp 3327–3337
#DOI: 10.1021/ac991500h

tof.h5 <- file.path('..','..','testdata/',"2015.07.17-10h40m34 Ethanol deurated Karl .h5")

fid <-H5Fopen(tof.h5)
tofblock <- get.raw.tofblock(fid)
indexhelp <- tof.indexhelp(tofblock)
cr <- make.curr.tofreader(tofblock, indexhelp)
#mc.table <- mass.calib.tof(tofblock, indexhelp)

kn <- 
  with(indexhelp, 
  	   if(N > 200)
			 floor(seq(from=1,to=N, length.out=100))
	  else	 seq(from=1,to=N)
)


locate.ions.block <- function(df.ion,indexhelp,
  preliminary.coeff=list(intercept = 900,square_mass = 17650 ),
  n.samp=100
  ){
  scl <- floor(seq(1,to=indexhelp$N, length.out=n.samp))
  pos.block <- expand.grid(scan=scl , ion=df.ion$ions)
  curry.read.line <- function (i) read.spec.ind(tofblock, indexhelp, i)
  pp <- apply(pos.block, 1,function (dd) { 
				curr.spec.line <- curry.read.line( dd[['scan']])
				pos <- find_ion_tof(dd[['ion']], preliminary.coeff, curr.spec.line)
  })

  pos.block$pos <- pp
  
}

# TODO: Massen suchen die *eindeutig* sind für die Kalibration

# h2oh3o = MonoisotopicMass(formula=list(O=1,x=1,H=4),charge=1, isotopes=list(x=17.999159)),
ions <- c(h3o = 21.0220875, 
		  #h2oh3o = MonoisotopicMass(formula=list(O=1,x=1,H=4),charge=1, isotopes=list(x=17.999159)),
		  h2o2h3o = MonoisotopicMass(formula=list(O=3,H=6),charge=1),
		  #no = MonoisotopicMass(formula=list(O=1,N=1))-electron.mass ,
		  #o2 = MonoisotopicMass(formula=list(O=2))-electron.mass ,
          aceton=59.04914,  
          isopren= MonoisotopicMass(formula=list(C=5,H=8),charge=1),
          etOH= MonoisotopicMass(formula=list(C=2,H=6,O=1),charge=1),
          butOH= MonoisotopicMass(formula=list(C=4,H=10,O=1),charge=1),
          butO= MonoisotopicMass(formula=list(C=4,H=8,O=1),charge=1),
          pentO= MonoisotopicMass(formula=list(C=5,H=10,O=1),charge=1),
          #hexO= MonoisotopicMass(formula=list(C=6,H=12,O=1),charge=1),
          ac= MonoisotopicMass(formula=list(C=2,H=4,O=1),charge=1)
          )
df.ion <- as.data.frame(ions)

pos.block <- locate.ions.block(df.ion,indexhelp, 
							   preliminary.coeff=list(intercept=1850, square_mass=17500))

pos.block <- pos.block[!is.na(pos.block$pos),]


pos.block <- within(pos.block,{
		 sq.mass <- sqrt(ion)
		 sq3.mass <- sq.mass*ion
		 sq5.mass <- sq3.mass*ion
  })
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
