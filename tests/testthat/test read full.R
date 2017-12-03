library(dplyr)
library(rhdf5)
library(devtools)

r.path <- file.path("..","..","R")

sl <- file.path(r.path,list.files(r.path, pattern="*.R"))
lapply(sl, source)

context("read full spectrum")

# cleanup -----------------------------------------------------------------

H5close()

base_path <-file.path("..","..","testdata")

# find files --------------------------------------------------------------

# list.files("testdata")
# tof.h5 <- file.path(base_path,"uncal.h5")
# h5ls(tof.h5)

# fm <- h5read(tof.h5,"FullSpectra/TofData")
# class(fm)


# h5ls(tof.h5)

test_that("open simple measurement successfully",{
  tof.h5 <- file.path(base_path,"Ac just exhale (2014-10-23T10h43m42_#).h5")
  fid <-H5Fopen(tof.h5)
  tofblock <- raw_tofblock(fid)
  indexhelp <- tof.indexhelp(tofblock)
  curr.spec.line <- read_spec_at(tofblock, indexhelp, 40)

  expect_equal(dim(curr.spec.line),398999, label= "length should always be 398999")
  expect_equal(curr.spec.line[1:10], as.array(rep(0,10)), label ="first 10 values are always 0")

  mc <- mass.calib.coeff.single(ions, preliminary.coeff=pars.approx, curr.spec.line)

  expect_equal(mc[["intercept"]], 900 , tolerance = 20, label = "intercept around 900")
  expect_equal(mc[["square_mass"]], 18000 , tolerance = 2000, label = "square coeff around 18000")

  # indexhelp$N
#   mass.calib <- data.frame(scan = (1:indexhelp$N) ) %>% rowwise() %>% do( {
#     curr.spec.line <- read_spec_at(tofblock, indexhelp, .$scan)
#     mc <- mass.calib.coeff.single(ions, preliminary.coeff=pars.approx, curr.spec.line)
#     data.frame(scan=.$scan, intercept=mc[["intercept"]], square_mass=mc[["square_mass"]])
#   } )

  H5close()
})

#
# plot(mass.calib$intercept, type="l")
# plot(mass.calib$square_mass, type="l")

# cc <- lm(cbind(intercept, square_mass) ~ ind, mc.fit)

# old ---------------------------------------------------------------------
#
#
# flat.tof <- read.tof.fullspectra(file.path("data","Ac just exhale (2014-10-23T10h43m42_#).h5"))
# mc.table <- mass.calib.tof(flat.tof)
# smooth.mc <- smooth.mass.cal(mc.table)
#
# mass.list <- c(21.022, 59.04914, 69.06988)
# mr1 <- lapply(mass.list, mass.to.range)
#
# mass.range <- mr1[[2]]
# c.ind <- 5
# curr.mass.cal <- smooth.mc[c.ind,]
#
# c.spec <- flat.tof[,c.ind]
# c.range <- mass.range2tof.range(curr.mass.cal, mass.range)
# ran <- c.range[["lower"]]:c.range[["upper"]]
# plot(c.spec[ran], type="l")
#
#
# one.peak.integrate(curr.mass.cal = c.mc, mass.range = cm, c.spec)
#
# pa <- apply(smooth.mc,1, function(x){
#   c.spec <- flat.tof[,x["ind"]]
#   one.peak.integrate(curr.mass.cal = x, mass.range = cm, c.spec)
# })
#
# plot(pa, type="l")
