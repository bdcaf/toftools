library(OrgMassSpecR)

nist.file = file.path('data', "Atomic Weights and Isotopic Compositions for All Elements.txt")

electron.mass <- 1/1820
isotopes <- list(C13= 13.0033548378, D=2.0141017780, N15=15.0001088984, O18=17.9991604)
mass <- MonoisotopicMass(formula = list(H=3, O=1), isotopes=list(O=isotopes$O18), charge=0) - electron.mass
mass

# H2O*H3o+
MonoisotopicMass(formula = list(H=5, O=1, x1=1), isotopes=list(x1=isotopes$O18)) - electron.mass

# 2*H2o*H3O+
library(OrgMassSpecR) - electron.mass

# 3*H2o*H3O+
MonoisotopicMass(formula = list(H=9, O=4)) - electron.mass

# NO+
MonoisotopicMass(formula = list(N=1, O=1)) - electron.mass

# Aceton+
MonoisotopicMass(formula = list(C=3, H=7, O=1)) - electron.mass

# Isopren+
MonoisotopicMass(formula = list(C=5, H=9)) - electron.mass


# Methacrolein+
MonoisotopicMass(formula = list(C=4, H=7, O=1)) - electron.mass
