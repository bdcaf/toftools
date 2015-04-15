context("read already integrated peaks")

tof.h5 <- file.path('..','..','data',"Ac just breath after C (2014-10-23T11h34m53_#).h5")
 
test_that("reading of the peak integrations",{
  m1 <- read.tof.peaks(tof.h5)
  expect_false(is.null(m1), info="not empty")
  expect_equal(c(3000,342), dim(m1),  info="correct size")
  expect_true( all(c('file',"(C5H5N)H+","41-(C3H5)+ (100.00%)") %in% colnames(m1)), info="check presence of some compounds")
})

