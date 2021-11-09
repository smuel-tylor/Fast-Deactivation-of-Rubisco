#A new version of Woodrow/Mott A*
#Using parameters for an FvCB A/Ci fit (Acifit12)
#to predict one-point Vcmax from A and Ci
#(after rearranging the FvCB eqn and substituting constants)

#requires the Rubisco.limited function from Acifit12
if (!exists("Rubisco.limited", mode = "function")) {
  source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/Acifit12.1.R"))
}

one.point.Vcmax <- function(Aci.mod, stoich, A, Pci, Pci.ss){
  
  Kco = Aci.mod$Kco
  Gamma.star = Aci.mod$Gamma.star
  Rd = Aci.mod$Rd
  gm = Aci.mod$gm
  
  Cc = Pci - A / gm
  
  Vc = (A + Rd) * (Cc + Kco) / (Cc - Gamma.star)
  
  Astar <- Rubisco.limited(Vcmax = Vc,
                           Kco,
                           Gamma.star,
                           Rd,
                           gm,
                           Pci.ss
                           )
  
  c(A, Pci, Vc, Astar)
  
}
