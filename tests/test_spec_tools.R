library(dplyr)
library(testthat)

source('R/spec_tools.R')
source('R/tofReader.R')
source('R/massCalTools.R')

fid <-H5Fopen(tof.h5)
tofblock <- get.raw.tofblock(fid)
indexhelp <- tof.indexhelp(tofblock)
aspec <- function(n) read.spec.ind(tofblock, indexhelp, n)

context("sparse spectrum")

sp0 <- aspec(1)

tb1 <- sparse2dense(sparse_spec(sp0))

expect_true(all(sp0 == tb1))
