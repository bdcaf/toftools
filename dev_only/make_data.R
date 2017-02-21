#' Example data set of ions used for mass calibration
#'
#' @name sample_ions
#' @docType data
#' @author Clemens Ager \email{clemens.ager@@uibk.ac.at}
#' @keywords data ions
NULL

electron.mass <- 1/1820

library(OrgMassSpecR)
ions <- c(h3o = 21.0220875, 
          h2oh3o = MonoisotopicMass(formula=list(O=1,x=1,H=4),charge=1, isotopes=list(x=17.999159)),
          h2o2h3o = MonoisotopicMass(formula=list(O=3,H=6),charge=1),
          no = MonoisotopicMass(formula=list(O=1,N=1))-electron.mass ,
          o2 = MonoisotopicMass(formula=list(O=2))-electron.mass ,
          acetone=59.04914,  
          isoprene= MonoisotopicMass(formula=list(C=5,H=8),charge=1),
          ethanol = MonoisotopicMass(formula=list(C=2,H=6,O=1),charge=1),
          butanol = MonoisotopicMass(formula=list(C=4,H=10,O=1),charge=1),
          butanone = MonoisotopicMass(formula=list(C=4,H=8,O=1),charge=1),
          pentanone = MonoisotopicMass(formula=list(C=5,H=10,O=1),charge=1),
          hexO= MonoisotopicMass(formula=list(C=6,H=12,O=1),charge=1),
          acectic_acid = MonoisotopicMass(formula=list(C=2,H=4,O=1),charge=1)
)
sample_ions <- as.data.frame(ions)

print(getwd())
save(sample_ions, file='data/sample_ions.Rdata')
