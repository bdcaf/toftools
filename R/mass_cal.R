#' initial fast mass calibration
#' @description uses the sum spectrum to estimate mass calibration
#' @param in_file file to quickly calibrate
#' @param masses list of masses to use for calibration
#' @param search_range in amu
#' @return lm fit of mass calibration
#' @export
mass_cal_fast <- function(in_file, masses, search_range = 0.3){
  sum_spec <- rawTof::get_sum_spec(in_file)
  mcfun <- read_stored_masscal(in_file)

  com_finder <- function(m_cur, wid = search_range) {
    range <- m_cur + c(-wid, wid)
    ind_range <- mcfun$to_index(range)
    pos_range <- seq(floor(ind_range[[1]]), ceiling(ind_range[[2]]))
    vals <- sum_spec[pos_range]
    com <- sum(pos_range * vals) / sum(vals)
    data.frame(mass = m_cur, index = com)
  }

  locations <- do.call(rbind, lapply(masses, com_finder))
  mass_calib <- lm(index ~ I(sqrt(mass)), data = locations)
}

#function(){
  #peak_list <- data_frame(center = unlist(masses),
                          #name = names(masses) ) %>%
  #mutate(
         #low = center - 0.1,
         #high = center + 0.1,
         #left_ind  = floor( predict(mc_fast, newdata = list(mass=low))),
         #right_ind = ceiling( predict(mc_fast, newdata = list(mass=high))),
         #wid_ind = right_ind - left_ind
         #)
  ##peak_list <- int_masses() %>% mutate(
                ##left_ind  = predict(mc_fast, newdata = list(mass=low)),
                ##right_ind = predict(mc_fast, newdata = list(mass=high)),
                ##wid_ind = right_ind - left_ind
                ##)

  #tof_ob <- new("TofClass", filename = in_file)
  #integration <- extract_tof(tof_ob, peak_list)
  #finalize(tof_ob)

  #with(dplyr::filter(integration, id == "h3o"),
       #plot(time, center_of_mass)
       #)
  #with(dplyr::filter(integration, id == "isoprene"),
       #plot(time, center_of_mass)
       #)
  #with(dplyr::filter(integration, id == "isoprene"),
       #plot(time, signal)
       #)
#}
