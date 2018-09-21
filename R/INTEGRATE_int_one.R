#' extract one measurement
#' @param tof_ob rawTof TofClass object
#' @param peak_list list as created by
#' @export
#' @import dplyr
#' @import rawTof
#'
#' @examples
#'
#' \dontrun{
#' tof_ob <- new("TofClass", filename = "raw/2018.04.20-09h17m18s_kuun.h5")
#' }

extract_tof <- function(tof_ob, peak_list) {
  timing <- rawTof::get_timing(tof_ob@.datafile)

  scan_row <- function(request) {
    w <- with(request, high - low) / 2
    centerd <- seq(-w, w, length.out = request$wid_ind)
    helpmat <- matrix(c(
      rep(1, length(centerd)),
      centerd,
      centerd^2,
      centerd^3
    ),
    byrow = TRUE, ncol = length(centerd)
    )
    block <- with(
      request,
      scan_block(tof_ob, ind = c(left_ind, wid_ind))
    )

    helpvals <- helpmat %*% block

    as.data.frame(t(helpvals)) %>%
      mutate(
        area = V1,
        e1 = V2 / area,
        e2 = V3 / area,
        e3 = V4 / area,
        center_of_mass = e1,
        variance = e2 - center_of_mass^2,
        skewness = e3 - 3 * e2 * e1 + 2 * e1^3,
        time = timing, id = request[["name"]]
      ) %>%
      select(id, time, signal = area, center_of_mass, variance, skewness)
  }

  peak_list %>%
    rowwise() %>%
    do(scan_row(.)) %>%
    ungroup() %>%
    mutate(numtime = as.double(time))
}


#' integrate one measurement from file
#' @param in_file hdf5 file
#' @param request list as created by
#' @export
#'
integrate_masses <- function(in_file, request)
  tof_wrap(in_file, extract_tof, request)
