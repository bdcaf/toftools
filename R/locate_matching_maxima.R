#' finds local maxima agreeing with masses, based on initial calibration
#' @param masses of calibration ions
#' @param sig signal to search in
#' @param mcfun previous mass calibration function
#' @note currently doesn't do any check whether a peak is present at the location. The user must take care that a peak is at the location.
#' @return tibble of local maxima
#' @import tibble
#' @export
locate_matching_maxima <- function(masses, sig, to_index){
  max_fun <- function(m_cur, wid = 0.06){
    range <- m_cur + c(-wid, wid)
    ind_range <- to_index(range)
    loc_sig <- sig[ ind_range[[1]]:ind_range[[2]]]
    c(mass = m_cur, index = which.max(loc_sig) - 1 + ind_range[[1]])
  }
  as_tibble(t(sapply(masses, max_fun)))
}
