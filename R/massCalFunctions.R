
# "i = p1*sqrt(m) + p2

.idx2mass <- function(idx, a, b){
  ((idx - b)  / a)^2
}



#' @title
#' locate ions for mass calibration
#' 
#' @description
#' searches for dominant ions in the relevant places and return a structure for mass calibration.
#' This should be improved soon to incorporate data alreade precalculated
#' 
#' 
#' @param df.ion named list of ion masses
#' @param scl scans to mass calibrate
#' @param perliminary.coef starting values for coefficients
#' @param n.samp number of points to take along the measurement
#' 
#' @return data.frame having columns scan, ion, and located position
locate.ions.block_ <- function(df.ion,
                               scl,
                               read.scan,
                               preliminary.coeff) {
  
  pos.block <- expand.grid(scan=scl , ion=df.ion$ions)
  pp <- apply(pos.block, 1,function (dd) { 
    curr.spec.line <- read.scan( dd[['scan']])
    pos <- find_ion_tof(dd[['ion']], preliminary.coeff, curr.spec.line)
  })
  
  pos.block$pos <- pp
  pos.block <- pos.block[!is.na(pos.block$pos),]  
  pos.block
}

#' @title
#' locate ions for mass calibration
#' 
#' @description
#' searches for dominant ions in the relevant places and return a structure for mass calibration.
#' This should be improved soon to incorporate data alreade precalculated
#' 
#' 
#' @param df.ion named list of ion masses
#' @param curr.h5 current tof measurement holder
#' @param perliminary.coef starting values for coefficients
#' @param n.samp number of points to take along the measurement
#' 
#' @return data.frame having columns scan, ion, and located position
#' 
#' @export
locate.ions.block <- function(df.ion,
                            curr.h5,
                            preliminary.coeff=list(intercept = 900,square_mass = 17650 ),
                            n.samp=100 ){
  with(curr.h5, {
    scl <- sampleScans(n.samp)
    locate.ions.block_(df.ion,scl,read.scan,preliminary.coeff)
  })
}

#' @title 
#' add variables to ions.block
#' 
#' @description
#' adds power 1/2, 3/2 and 5/2 of the masses to the position block in order for fitting the more complex mass calibration function of Muller 2013 [1]
#' @references
#' [1] M. Müller, T. Mikoviny, W. Jud, B. D’Anna, A. Wisthaler, and B. D’Anna, “A New Software Tool for the Analysis of High Resolution PTR-TOF Mass Spectra,” Chemom. Intell. Lab. Syst., vol. 127, pp. 158–165, Jun. 2013.
#' @export
enrich.ions.block <- function(pos.block){
  within(pos.block,{
    sq.mass <- sqrt(ion)
    sq3.mass <- sq.mass*ion
    sq5.mass <- sq3.mass*ion
  })  
}

#' @title 
#' fit mass calibration function
#' 
#' @description
#' uses the more elaborate function described in Muller 2013 [1]
#' 
#' @param curr.h5 current measurement
#' @param pos.block enriched position data of ions
#' @param num.knots number of knots for spline
#' 
#' @return linear model fit
#' 
#' @references
#' [1] M. Müller, T. Mikoviny, W. Jud, B. D’Anna, A. Wisthaler, and B. D’Anna, “A New Software Tool for the Analysis of High Resolution PTR-TOF Mass Spectra,” Chemom. Intell. Lab. Syst., vol. 127, pp. 158–165, Jun. 2013.
#' @export
fit.mass.calib <- function(curr.h5,pos.block,num.knots=3){
  kn <- floor(seq(1,to=curr.h5$num.scan(), length.out=num.knots+2))[2:num.knots+1]  
  rlm(pos ~ ion + sq3.mass + sq5.mass + sq.mass * ns(scan,knots=kn), pos.block, method='MM')
}


#' prepare for index calculations of block
#' 
#' given a required list of ion masses the relevant indices at given scan is calculated
#' @param curr.h5 the h5 for which the calculation is to be performed
#' @param masslist masses to be considered
#' @return data.frame of prepared variables
#' @export
#' @examples
#' \dontrun{
#' masslist <- seq(46.9,47.1, by=0.001)
#' make.read.block(curr.h5,masslist)
#' }
make.read.block <- function(curr.h5,masslist, scans=seq(curr.h5$num.scan())){  
  dfm <- expand.grid(ion=masslist, scan=scans)
  enrich.ions.block(dfm)
}


#' predict position index of masses
#' 
#' uses the parameters of readblock to predict index of specified masses
#' @export
#' @param fitted model
#' @param required masses
#' @return data.frame of index position
#' @examples
#' \dontrun{
#' masslist <- seq(46.9,47.1, by=0.001)
#' rmod <- fit.mass.calib(curr.h5, pos.block)
#' rb <- make.read.block(curr.h5,masslist)
#' pred.pos <- predic.positions(rmod, rb)
#' }
predic.positions <- function(model,read.block){
  pm <- predict(model, newdata=read.block)
  with(read.block, 
       data.frame(scan=scan, ion=ion, index=pm))
}