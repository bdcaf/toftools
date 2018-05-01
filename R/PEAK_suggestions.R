
rle_clean_step <- function(rl, min_l, dir){
  bad <- with(rl, lengths < min_l & values == dir)
  rl$values[bad] <- !rl$values[bad]
  ig2 <- inverse.rle(rl)
  rle(ig2)
}

#' no longer used
rle_clean_rec <- function(rl, min_l, dir = TRUE){
  old <- NA
  while (old != length(rl$values)){
    old <- length(rl$values)
    rl <- rle_clean_step(rl, min_l, dir)
  }
  rl
}

#' clean rle
#' @import magrittr
rle_cleaner <- function(rl, min_l){
  rle_clean_step(rl, min_l, dir = T) %>%
    rle_clean_step(min_l, dir = F)
}

#' calculate center of mass
com <- function(ss, left, right){
  ix <- left:right
  cs <- ss[ix]
  nc <- sum(cs)
  mc <- sum(cs * ix)
  mc/nc
}

#' suggest peaks base on sum spectrum and mass calib
#' @param sumspec sum spectrum
#' @param masscal mass calibration from lm
#' @param min_sig minimal signal in sum spectrum to consider a peak
#' @param min_wid minimal width of a peak
#' @description this currently doesn't do any deconvolution of peaks
#' @import dplyr
#' @export
peak_suggestions <- function(sumspec, masscal,
                             start_mass = 20.5,
                             min_sig = 1e4,
                             min_wid = 5){
  ind_start <- ceiling(predict(masscal,
                               newdata = list(mass = start_mass)))
  ig <- sumspec > min_sig & seq_along(sumspec) > ind_start
  s_rl <- rle(as.vector(ig)) %>% rle_cleaner(min_l = min_wid)


  inv_mc <- function(ind) {
    ((ind - masscal$coefficients[[1]])/(masscal$coefficients[[2]]))^2
  }

  data_frame(lengths = s_rl$lengths, values = s_rl$values) %>%
    mutate(right_ind = cumsum(lengths),
           left_ind = lag(right_ind),
           wid_ind = right_ind - left_ind) %>%
  dplyr::filter(values) %>%
  rowwise() %>%
  mutate(center_ind = com(ss, left_ind, right_ind)) %>%
  ungroup() %>%
  mutate(center = inv_mc(center_ind),
         name = sprintf("m/z %.2f", center)
         )
}
