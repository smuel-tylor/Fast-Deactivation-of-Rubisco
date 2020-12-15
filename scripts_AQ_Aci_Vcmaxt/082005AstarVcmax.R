#A new version of Woodrow / Mott A*
#here using A and ci to predict Vcmax for each single point by rearranging the FvCB eqn and substituting in constants

#requires the Rubisco.limited function from
source("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rscripts/081913Acifit12.1.R")

one.point.Vcmax<-function(Aci.mod,stoich,A,Pci,Pci.ss){

Kco=Aci.mod$Kco
Gamma.star=Aci.mod$Gamma.star
Rd=Aci.mod$Rd
gm=Aci.mod$gm

Cc=Pci-A/gm

Vc=(A+Rd)*(Cc+Kco)/(Cc-Gamma.star)

Astar<-Rubisco.limited(Vcmax=Vc,Kco,Gamma.star,Rd,gm,Pci.ss)

c(A,Pci,Vc,Astar)

}