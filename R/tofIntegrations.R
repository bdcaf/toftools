#' helper function to create integration borders as in TofDAQ Viewer
#' @export
mass.to.range <- function(mass=21.022, res=1000) c(center=mass, lower=mass-mass/res, upper=mass+mass/res)

#' helper to convert a mass range into tof timing
mass.range2tof.range <- function(curr.mass.cal, mass.range){
  left <- floor(   curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']]*sqrt(mass.range[['lower']]))
  right <- ceiling(curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']]*sqrt(mass.range[['upper']]))
  c(lower=left, upper=right)
}

#' integrate single peak in single spec
one.peak.integrate <- function (curr.mass.cal, mass.range, c.spec) {  
  c.range <- mass.range2tof.range(curr.mass.cal, mass.range)
  ran <- c.range[['lower']]:c.range[['upper']]
  sum(c.spec[ran])
}

#' integrate one compound over tof dataset
#' @export
integrate.peak.over.set <- function(mass.cal.list, mass.range, flat.tof){
  apply(mass.cal.list,1, function(x){  
    c.spec <- flat.tof[,x['ind']]
    one.peak.integrate(curr.mass.cal = x, mass.range = mass.range, c.spec)
  })
}