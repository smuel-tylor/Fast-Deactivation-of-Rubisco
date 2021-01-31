#Master script for Gonzalez Escobar et al.
#runs all analysis scripts in correct sequence

library(here)

source(here("1_scripts_response_curves/082005summarize.R"))

source(here("2_scripts_sun_shade_sun/082005NaturePlantsnlmeVcmax.R"))

source(here("2_scripts_sun_shade_sun/082005NaturePlantsnlmeLightResponses.R"))

source(here("2_scripts_sun_shade_sun/082005NaturePlantsnlmeActivationState.R"))

source(here("3_scripts_diurnal_cowpea/082005NaturePlantsdiurnalmodels.R"))

#Figures
source(here("2_scripts_sun_shade_sun/082005NaturePlantsFigure1.R"))

source(here("3_scripts_diurnal_cowpea/082005NaturePlantsFigure2.R"))
