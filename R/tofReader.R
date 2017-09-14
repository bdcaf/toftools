library(rhdf5)
library(dplyr)
library(tidyr)
library(purrr)
library(Matrix)

#' reads integrated peaks as integrated by PTR-MS TOD-DAQ-Viewer.
#'
#' WARNING: these must be first integrated in the software otherwise and error
#' will be thrown.
#' @export
#' @param tof.h5 filename of hdf5 file
#' @return dataframe with structure: rows = scans, cols = counts + information
#' @examples
#' tof.h5 <- file.path("data","uncal.h5")
read.tof.peaks <- function( tof.h5 ){
  fid <-H5Fopen(tof.h5)
  at <- H5Aopen(fid, "NbrWaveforms")
  waveforms <- H5Aread(at)
  H5Aclose(at)
  pd1 <- H5Dopen(fid, "PeakData/PeakTable")
  peak.data <- H5Dread(pd1)
  pd2 <- H5Dopen(fid, "PeakData/PeakData") # peak value in counts per extractions
  peak.value <- H5Dread(pd2)

  peak.counts <- peak.value * waveforms[[1]]
  tmp <- apply(peak.counts, 1, as.vector)
  peak.frame <- as.data.frame(tmp)
  colnames(peak.frame) <- peak.data$label
  peak.frame$file <- tof.h5
  return(peak.frame)
}

read.sum.spec <- function(fid){
  gr <- H5Dopen(fid, "FullSpectra/SumSpectrum")
  H5Dread(gr)
}
read.mass.cals <- function(fid){
  gr <- H5Dopen(fid, "FullSpectra/MassCalibration")
  H5Dread(gr)
}

#' approximate mass calibration
pars.approx <- list(intercept = 900,square_mass = 17600 )

#' sample ions for TOF calibration
ions <- c(h3o = 21.0220875,
          no=29.99744,
          o2 = 31.989281,
          aceton=59.04914,
          h5o2=39.03265,
          h7o3=55.03897,
          meth_formate = 61.028406,
          methacrolein = 71.04914)

#' calculate mass from TOF index
mass.calc <- function(pars, vec) {(( vec - pars[[1]] )/ pars[[2]] )^2 }

#' find location of ion maximum in TOF spectrum
find_ion_tof <- function (ion.mass = 21.0220875, preliminary.coeff, curr.spec.line, wid=0.3) {
  mass.range <- round(preliminary.coeff$intercept + preliminary.coeff$square_mass*sqrt(ion.mass+ wid*c(-1,1)))
  ig <- mass.range[1]:mass.range[2]
  # plot(curr.spec.line[ig], type="l")
  pos <- ig[1]-1 + which.max(curr.spec.line[ig])
  ifelse(curr.spec.line[pos] > 10 & curr.spec.line[pos] < 1e6, pos, NA) # NA if either saturated or not enough
}

mass.calib.coeff.single <- function (ions, preliminary.coeff, curr.spec.line) {
  tofs <- sapply(ions, function(x) find_ion_tof(x, preliminary.coeff, curr.spec.line))
  mf <- data.frame(sqmass=sqrt(ions), tofs = (tofs))
  mod <- lm(tofs ~  sqmass, mf)
  co <- coef(mod)
  out <- data.frame(intercept = co[[1]], square_mass = co[[2]])
  return(out)
}

#' calculates mass calibration coefficients for every scan in the flat tof array
#' @note this thkes about 30s for a measurement containing 3000 scans
#' @param tofblock from get.raw.tofblock
#' @param indexhelp from tof.indexhelp
#' @param preliminary.coeff some start values for coefficients
#' @param n.ions ion list to be used for calibration
#' @return data frame of calibration coefficients
#' @examples
#' tof.h5 <- file.path("data","Ac just breath after C (2014-10-23T11h34m53_#).h5")
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' cr <- make.curr.tofreader(tofblock, indexhelp)
#' mc.table <- mass.calib.tof(tofblock, indexhelp)
#' @export
mass.calib.tof <- function(tofblock, indexhelp,
                           preliminary.coeff=list(intercept = 900,square_mass = 17650 ),
                           n.ions = unlist(ions)){
  data.frame(scan = (1:indexhelp$N) ) %>% rowwise() %>% do( {
    curr.spec.line <- read.spec.ind(tofblock, indexhelp, .$scan)
    mc <- mass.calib.coeff.single(ions, preliminary.coeff=pars.approx, curr.spec.line)
    data.frame(scan=.$scan, intercept=mc[["intercept"]], square_mass=mc[["square_mass"]])
  } )
}

#' smoothes the mass calibration to avoid local jumps
#' @param mc.table raw mass calibration
#' @return mass calibration table after smoothing
#' @examples
#' mc.table <- mass.calib.tof(tofblock, indexhelp)
#' smooth.mc <- smooth.mass.cal(mc.table)
smooth.mass.cal <- function (mc.table) {
  mc.fit <- mc.table
  mc.fit$ind <- 1:nrow(mc.fit)
  cc <- lm(cbind(intercept, square_mass) ~ ind, mc.fit)
  pv <- predict(cc,  data.frame(ind=1:nrow(mc.fit)))
}

#' first step of reading - get the tof block
#' not export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' H5close(fid)
get.raw.tofblock <- function(fid) {
  gr <- H5Gopen(fid,"FullSpectra")
  H5Dopen(gr,"TofData")}

#' returns full block of tof spectra
#'
#' @descriptions reads the full block of tof spectra
#' This can be very large!
#' @param tofblock from get.raw.tofblock
#' @return array, first dimenstion is timebin, second scan index
get.full.specblock <- function(tofblock){
  h5spaceFile <- H5Dget_space(tofblock)
  dims <- H5Sget_simple_extent_dims(h5spaceFile)
  ot <- H5Dread(h5dataset = tofblock,
          h5spaceFile = h5spaceFile,
          compoundAsDataFrame = FALSE)
  #ot2 <- ot
  #dd <- dim(ot)
  #dim(ot2) <- c(dd[[1]], prod(dd[-1]))
  #all(ot2[,1]==ot[,1,1,1])
  #all(ot2[,2]==ot[,1,2,1])
  dd <- dim(ot)
  dim(ot) <- c(dd[[1]], prod(dd[-1]))
  Matrix(ot, sparse=T)
}

#' pre calculates some values
#' not export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
tof.indexhelp <- function(tofblock){
  h5spaceFile <- H5Dget_space(tofblock)
  dims <- H5Sget_simple_extent_dims(h5spaceFile)
  calc.indices <- expand.grid(buf=1:dims$size[[3]], write=1:dims$size[[4]])
  h5spaceMem <- H5Screate_simple(dims$size[[1]])
  list(N=nrow(calc.indices),
       h5spaceFile=h5spaceFile,
       h5spaceMem=h5spaceMem,
       calc.indices=calc.indices,
       dims=dims$size)
}

#' reads the spec of a single tof scan
#' @export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' spec <- read.spec.ind(tofblock, indexhelp, 40)
#' plot(spec,type="l")
read.spec.ind <- function(tofblock, indexhelp, i){
  pos <- indexhelp$calc.indices[i,]
  H5Sselect_hyperslab(indexhelp$h5spaceFile,
                      start = c(1,1,pos$buf,pos$write),
                      count = c(indexhelp$dims[[1]],1,1,1) )
  H5Dread(h5dataset = tofblock,
          h5spaceFile = indexhelp$h5spaceFile,
          h5spaceMem = indexhelp$h5spaceMem )
}

#' helper to wrap reader call
#' not export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' cr <- make.curr.tofreader(tofblock, indexhelp)
#' spec5 <- cr(5)
make.curr.tofreader <- function(tofblock, indexhelp){
  function(i){read.spec.ind(tofblock,indexhelp,i)}
}

h5ToSqlite <- function(path, tof.h5){
  src <- src_sqlite(path, create=T)
  #sparse.table <- tbl("sparse", src)

  fid <-H5Fopen(tof.h5)
  tofblock <- get.raw.tofblock(fid)
  indexhelp <- tof.indexhelp(tofblock)
  aspec <- function(n) read.spec.ind(tofblock, indexhelp, n)

  #db_drop_table(src$con, "sparse")

  ind <- 4

  osframe <- function(ind) {
    as <- aspec(ind)
    ss <- sparse_spec(as)
    cbind(scan = ind, ss)
  }

    #ss2 <- ss %>%
      #rowwise() %>%
      #mutate(ser = list(serialize(v,NULL))) %>%
      #select(-v)
  system.time(
  tmp <- data_frame(scan = seq_len(indexhelp$N)) %>%
    rowwise() %>%
    do( osframe(.$scan))
    )

  to_db <- tmp %>%
      rowwise() %>%
      mutate(ser = list(serialize(v,NULL))) %>%
      select(-v)

  copy_to(src, to_db, name="sparse", temporary=F)
  #copy_to(src, tmp, name="sparse", temporary=F, indexes=list("scan","bin"))
}


###############################
# S3 class for easier working #
###############################

#' tof_h5 constructor for S3 class
#' @export
tof_h5 <- function(tof.h5){
  fid <-H5Fopen(tof.h5)
  tofblock <- get.raw.tofblock(fid)
  indexhelp <- tof.indexhelp(tofblock)
  structure( list(filename = tof.h5,
                  fid = fid,
                  tofblock = tofblock,
                  indexhelp = indexhelp),
            class = "tof_h5")
}

#' print tof_h5 class
#' @export
print.tof_h5 <- function(th5){
  with(th5, cat(sep="\n",
                "H5 TOF measurement",
                paste("file:",filename),
                paste("scans:",indexhelp$N)
                ))
}


#' reads the spec of a single tof scan from tof_h5
#' @export
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' spec <- read.spec.ind(tofblock, indexhelp, 40)
#' plot(spec,type="l")
read_spec_ind.tof_h5 <- function(tofH, i){
  with(tofH,{
  pos <- indexhelp$calc.indices[i,]
  H5Sselect_hyperslab(indexhelp$h5spaceFile,
                      start = c(1,1,pos$buf,pos$write),
                      count = c(indexhelp$dims[[1]],1,1,1) )
  H5Dread(h5dataset = tofblock,
          h5spaceFile = indexhelp$h5spaceFile,
          h5spaceMem = indexhelp$h5spaceMem ) }
  )}


sum_spec.tof_h5 <- function(tofH)
  with(tofH, {
         gr <- H5Dopen(fid, "FullSpectra/SumSpectrum")
         H5Dread(gr)
  })

stored_mass_cal.tof_h5 <- function(tofH){
  gr <- h5readAttributes(tofH$fid, "FullSpectra")
  with(gr,
       list( to_mass = Vectorize(function(i) ((i-`MassCalibration p2`)/`MassCalibration p1`)^2),
             to_index = Vectorize(function(m) `MassCalibration p2` + `MassCalibration p1`*sqrt(m))
             ))
}
