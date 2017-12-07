#' helper function to create a reusable function for warping
#' @export
#' @param v spectrum
#' @param ind position of v
#' @return warp function that recalculates the spectrum
make_densewarp <- function(v,  inds = seq_along(v)){
  vc <- local({
    cuv <- cumsum(v)
    lastcu <- last(cuv)
    c(0,  cuv,  lastcu)
  })
  index_along <- c(inds[[1]] - 1, inds, inds[[length(inds)]] + 1)
  function(warp_fun, where = inds){
    i_trans <- warp_fun(index_along)
    i_new <- c(where[[1]] - 1, where)
    v_new <- diff( approx(i_trans,  vc,  i_new,
                          yleft = 0,
                          yright = last(vc))$y)
    return(v_new)
  }
}
