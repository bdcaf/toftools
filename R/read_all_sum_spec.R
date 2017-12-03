#' read all sum spec s
#' @import magrittr
#' @importFrom rawTof get_sum_spec
#' @param data_files list of h5 files
#' @return list of sum spectrum vectors
#' @export
read_all_sum_spec <- function(data_files){
  init_spec <- lapply(data_files, get_sum_spec)
  len <- sapply(init_spec, length) %>% min
  sum_specs <- lapply(init_spec, function(x) x[1:len])
}
