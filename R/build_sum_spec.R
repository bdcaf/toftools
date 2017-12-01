#' average sum spec
#' @param data_files list of h5 files
#' @return sum spectrum vector
#' @export
build_sum_spec <- function(data_files){
  sum_specs <- read_all_sum_spec(data_files)
  rr <- do.call(rbind, sum_specs)
  apply(rr, 2, mean)
}
