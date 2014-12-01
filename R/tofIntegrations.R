library(dplyr)

#' adds mass bounds to mass list
#' @export
#' @examples
#' mass.list <- read.table(text='
#'  compound  center
#'  H3O+  21.022
#'  acetone   59.049
#'  isoprene  69.0699',
#' stringsAsFactor=FALSE, header=TRUE)
#' create.mass.bounds(mass.list)
create.mass.bounds <- function(mass.list, res=1500){
  mass.list %>% mutate(lower = center-center/res, upper = center+center/res)
}

#' adds mass indices to mass list
#' @export
#' @examples
#' mass.list <- read.table(text='
#'  compound  center
#'  H3O+  21.022
#'  acetone   59.049
#'  isoprene  69.0699',
#' stringsAsFactor=FALSE, header=TRUE)
#' mass.list <- create.mass.bounds(mass.list)
#' add.index.range(mass.list, curr.mass.cal)
add.index.range <- function(mass.list, curr.mass.cal){
  mutate(mass.list, 
         lower.index = curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']] * sqrt(lower),
         upper.index = curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']] * sqrt(upper)) %>%
    mutate(upper.index = ceiling(upper.index), lower.index = floor(lower.index))
}

#' integrate single spec line
#' @examples
#' mass.list <- read.table(text='
#'  compound  center
#'  H3O+  21.022
#'  acetone   59.049
#'  isoprene  69.0699',
#' stringsAsFactor=FALSE, header=TRUE)
#' mass.list <- create.mass.bounds(mass.list)
#' mass.list <- add.index.range(mass.list, curr.mass.cal)
#' integrate.one.spec(mass.list2, spec)
integrate.one.spec <- function(mass.list, spec){
  rowwise(mass.list) %>% mutate(area = sum(spec[lower.index:upper.index]))%>%
    select(compound, center, area)
}

#' @export
#' integrates the full spec
#' @examples
#' fid <-H5Fopen(tof.h5)
#' tofblock <- get.raw.tofblock(fid)
#' indexhelp <- tof.indexhelp(tofblock)
#' curr.reader <- make.curr.tofreader(tofblock, indexhelp)
#' mass.list <- read.table(text='
#'  compound  center
#'  H3O+  21.022
#'  acetone   59.049
#'  isoprene  69.0699',
#' stringsAsFactor=FALSE, header=TRUE)
#' mass.list <- create.mass.bounds(mass.list)
#' int.res <- integrate.full.spec(curr.reader, indexhelp, smooth.mc, mass.list)
#' library(ggplot2)
#' ggplot(int.res) + geom_line(mapping=aes(x=index, y=area, color=compound)) + scale_y_log10()
integrate.full.spec <- function(curr.reader, indexhelp, mass.cal, mass.list){
  target <- data.frame(mass.cal) %>% mutate(index = 1:indexhelp$N)
  
  rowwise(target) %>% do({
    mass.list2 <- add.index.range(mass.list2, .)
    spec <- curr.reader(.$index)
    res <- integreate.one.spec(mass.list2, spec)
    res$index=.$index
    return(res)
  })
}