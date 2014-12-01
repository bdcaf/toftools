fid <-H5Fopen(tof.h5)
tofblock <- get.raw.tofblock(fid)
indexhelp <- tof.indexhelp(tofblock)
cr <- make.curr.tofreader(tofblock, indexhelp)
mc.table <- mass.calib.tof(tofblock, indexhelp)
smooth.mc <- smooth.mass.cal(mc.table)


# integrate one peak ------------------------------------------------------
mass.list <- read.table(text='
                        compound  center
                        H3O+  21.022
                        acetone   59.049
                        isoprene  69.0699',
                        stringsAsFactor=FALSE, header=TRUE)

res = 1500
mass.list %>% mutate(lower = center-center/res, upper = center+center/res)


mass.list2 <- mass.list %>% create.mass.bounds(res=1000)

mass.list2 <- mass.list2 %>% mutate(lower.index = curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']] * sqrt(lower),
                      upper.index = curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']] * sqrt(upper)) %>%
  mutate(upper.index = ceiling(upper.index), lower.index = floor(lower.index))

mass.list2 <- add.index.range(mass.list2, curr.mass.cal)

mass.list2 %>% rowwise() %>% mutate(area = sum(spec[lower.index:upper.index])) %>%
  select(compound, center, area)

integreate.one.spec(mass.list2, spec)


# over whole sample -------------------------------------------------------

mass.list <- read.table(text='
                        compound  center
                        H3O+  21.022
                        acetone   59.049
                        isoprene  69.0699',
                        stringsAsFactor=FALSE, header=TRUE)
mass.list <- mass.list %>% create.mass.bounds(res=1000)

target <- data.frame(smooth.mc) %>% mutate(index = 1:indexhelp$N)

t2 <- rowwise(target) %>% do({
  mass.list2 <- add.index.range(mass.list2, curr.mass.cal)
  spec <- cr(.$index)
  res <- integreate.one.spec(mass.list2, spec)
  res$index=.$index
  return(res)
  })


library(ggplot2)
ggplot(t2) + geom_line(mapping=aes(x=index, y=area, color=compound)) + scale_y_log10()


data.frame(index = 1:indexhelp$N) %>% rowwise() %>% do( {} )


c.ind <- 5
curr.mass.cal <- smooth.mc[c.ind,]
spec <- cr(c.ind)

mass.range <- mr1[[2]]
curr.index.range <- curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']] * sqrt(mass.range)
plot(spec[ curr.index.range[['lower']]:curr.index.range[['upper']] ], type='l')


area <- sum(spec[ curr.index.range[['lower']]:curr.index.range[['upper']] ])

one.peak.integrate(curr.mass.cal, mass.range, spec)
 

# several peaks -----------------------------------------------------------
# one.spec.full.int <- function (mr1, curr.mass.cal, spec) {
#   lapply(mr1, function(cm){
#     curr.index.range <- curr.mass.cal[['intercept']] + curr.mass.cal[['square_mass']] * sqrt(cm)
#     sum(spec[ curr.index.range[['lower']]:curr.index.range[['upper']] ])  
#   })
# }

one.spec.full.int(mr1, curr.mass.cal, spec)
