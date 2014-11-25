library(dplyr)
library(tidyr)

tof.h5 <- file.path('data',"uncal.h5")
h5ls(tof.h5)

fm <- h5read(tof.h5,'FullSpectra/TofData')
class(fm)

preliminary.coeff <- list(intercept = 900,square_mass = 17650 )
pars <- pars.approx

# plot(ma, type='l')

aa <- array(fm, dim=c(dim(fm)[1],prod(dim(fm)[2:4]))) # reshape
aa1 <- data.frame(aa)
colnames(aa1) <- 1:120
a1 <- as.tbl(aa1) %>% mutate(.,tof = as.numeric(row.names(.))) %>% gather(scan,val,-tof) %>% mutate(scan=as.numeric(scan))

n.ions <- unlist(ions)
mass_calibs <- a1 %>% group_by(scan) %>% do( {
  mass.calib.coeff.single(n.ions = n.ions, preliminary.coeff=preliminary.coeff, curr.spec.line = .$val)
} ) # takes some time

# ap.a <- apply(aa,2, function(x) print(length(x)))

ap.a <- apply(aa,2, function(x) mass.calib.coeff.single(n.ions = n.ions, preliminary.coeff=preliminary.coeff, curr.spec.line = x))
df.a <- as.data.frame(matrix(unlist(ap.a), ncol=2, byrow=TRUE))
colnames(df.a) <- c('intercept','square_mass')
df.a


curr.spec.line <- aa[,1]
mass.calib.coeff.single(n.ions = n.ions, preliminary.coeff=preliminary.coeff, curr.spec.line = curr.spec.line)




a1 <- as.tbl(aa)

cs <- aa[,1]
ma <- mass.calc (pars.approx,1:length(cs))
plot(ma,cs, type='l')

# aa %>% gather(tof,-X1)

tmp <- apply(fm, 1, as.vector) 


# plot(ma, tmp[3,], type='l')



ci <- ions[[1]]
cs <- aa[,1]
plot(preliminary.mass[ig], cs[ig], type='l')


find_ion_tof <- function (ion.mass = 21.0220875, preliminary.mass, spectrum) {  
  ig <- preliminary.mass > ion.mass-0.3 & preliminary.mass < ion.mass+0.3  
  pos <- match(TRUE,ig) + which.max(spectrum[ig])
  ifelse(cs[pos] > 30, pos, NA)
}

find_ion_tof(ions[[2]], ma, cs)
plot(aa[97563,], type='l')
# 
# ions <- list(h3o = 21.0220875, no=29.99744, aceton=59.04914, isopren=69.06988, h5o2=39.03265, h7o3=55.03897)
n.ions <- unlist(ions)
mcts <- mass.calib.tof.single(n.ions = n.ions, preliminary.mass=ma, curr.spec.line = cs)
plot(mcts, type='l')

apply(tmp,1, length)

 
out <- data.frame(mass = mass.calc (co,1:ncol(tmp)), signal=cs)
ig <- out$mass >68 & out$mass<70
plot(out[ig,], type='l')
