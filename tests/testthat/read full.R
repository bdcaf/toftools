tof.h5 <- file.path('data',"uncal.h5")
h5ls(tof.h5)

fm <- h5read(tof.h5,'FullSpectra/TofData')

pars.approx <- list(b = 900,a = 17600 )
pars <- pars.approx

mass.calc <- function(pars, vec) {(( vec - pars[[1]] )/ pars[[2]] )^2 }

ma <- mass.calc (pars.approx,1:ncol(tmp))

plot(ma, type='l')

tmp <- apply(fm, 1, as.vector) 


plot(ma, tmp[3,], type='l')



ci <- ions[[1]]
cs <- tmp[3,]
plot(preliminary.mass[ig], cs[ig], type='l')

find_ion_tof <- function (ion.mass = 21.0220875, preliminary.mass, spectrum) {  
  ig <- preliminary.mass > ion.mass-0.3 & preliminary.mass < ion.mass+0.3  
  pos <- match(TRUE,ig) + which.max(spectrum[ig])
  ifelse(cs[pos] > 30, pos, NA)
}

# find_ion_tof(ions[[2]], ma, cs)

ions <- list(h3o = 21.0220875, no=29.99744, aceton=59.04914, isopren=69.06988, h5o2=39.03265, h7o3=55.03897)
n.ions <- unlist(ions)
tofs <- sapply(n.ions, function(x) find_ion_tof(x, ma, cs))

plot(ions, log10( tofs))

pos <- tofs[[6]]
plot(cs[(pos-100):(pos+100)], type='l')


mf$tofs2 <- mf$tofs^2
mf$sqmass <- sqrt(mf$mass)
mf
mod <- lm(mass ~ tofs2 + tofs, mf)
mod
plot(mod)

anova(mod)

mf <- data.frame(sqmass=sqrt(n.ions), tofs = (tofs))
mod <- lm(tofs ~  sqmass, mf)
co <- coef(mod)

out <- data.frame(mass = mass.calc (co,1:ncol(tmp)), signal=cs)

plot(mod)
anova(mod)
