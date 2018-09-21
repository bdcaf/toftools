rhdf5::H5close()
rhdf5::h5closeAll()
context("read already integrated peaks")

tof.h5 <- file.path("../..","testdata/2015.07.17-10h40m34 Ethanol deurated Karl .h5")

test_that("reading of the peak integrations",{
  expect_true(file.exists(tof.h5),
              info = "check test file presence")

  m1 <- read.tof.peaks(tof.h5)

  expect_equal(c(141,499), dim(m1),  info="correct size")
  expect_true( 'file' %in% colnames(m1),
              info = "fileinfo present.")
  expect_true( 'file' %in% colnames(m1),
              info = "fileinfo present.")
  expect_false( any(duplicated(colnames(m1))),
              info = "test for duplicated column names")
})
