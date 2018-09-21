#' wraps opening and closing of tof object
#' @param in_file hdf5 file
#' @param to_function tofob function that takes object in first plase
#' @param ... furthetr param
#' @import rawTof
#' @export
#'
tof_wrap <- function(in_file, toFunction, ...){
    tof_ob <- new("TofClass", filename = in_file)
    on.exit( finalize(tof_ob) )
    toFunction(tof_ob, ...)
  }
