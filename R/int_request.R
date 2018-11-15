#' Function to create list of locations in tof file
#' @import dplyr
#' @export
#' @param masses list of masses to integrate
#' @param fitm mass calibration that can predict the location in original frame
#' @return tibble with locations of peaks in tof frame
request_list <- function(masses, fitm)
  tibble(name = names(masses), mass = unlist(masses)) %>%
    mutate(
      wid = sqrt(mass) / 300, low = mass - wid, high = mass + wid,
      left_ind = predict(fitm, newdata = data.frame(mass = low)) %>% floor(),
      right_ind = predict(fitm, newdata = data.frame(mass = high))
      %>%
        ceiling(),
      wid_ind = right_ind - left_ind
    )

#' Function to create list of locations in tof file
#' @import dplyr
#' @export
#' @param masses list of masses to integrate
#' @param function to that can predict the location in original frame
#' @return tibble with locations of peaks in tof frame
request_mass <- function(masses, mass_fun)
  tibble(name = names(masses), mass = unlist(masses)) %>%
    mutate(
      wid = sqrt(mass) / 300, low = mass - wid, high = mass + wid,
      left_ind = mass_fun(low) %>% floor(),
      right_ind = mass_fun(high) %>% ceiling(),
      wid_ind = right_ind - left_ind
    )
