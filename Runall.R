#Master script for Gonzalez Escobar et al.
#runs all analysis scripts in correct sequence
#see readme.txt for further details

library(here)

#note this script calls cleanup.R
source(here("1_scripts_response_curves/summarize.R"))

source(here("2_scripts_sun_shade_sun/nlmeVcmax.R"))

source(here("2_scripts_sun_shade_sun/nlmeLightResponses.R"))

source(here("2_scripts_sun_shade_sun/nlmeActivationState.R"))

source(here("3_scripts_diurnal_cowpea/diurnalmodels.R"))

#Figures
source(here("2_scripts_sun_shade_sun/Figure1_Multipanel_S_Vcmax_Vi_Vt.R"))

source(here("3_scripts_diurnal_cowpea/Figure2_DiurnalModel.R"))

#Supplementary Figures
#Figure S1 was a schematic, not associated with data analysis

#Figures S2-S4 were composites of the graphs generated below + schematics
source(here("5_scripts_light_rig_supplements/FigureS2_DirectVsIncubation.R"))

source(here("5_scripts_light_rig_supplements/FigureS3_IncubationTime.R"))

source(here("5_scripts_light_rig_supplements/FigureS4_Bicarbonate.R"))

source(here("5_scripts_light_rig_supplements/FigureS5_LightCurves.R"))

#Figure S6 was a photograph

source(here("1_scripts_response_curves/FigureS7_AciResponses.R"))

source(here("2_scripts_sun_shade_sun/FigureS8_LimitationStates.R"))

source(here("2_scripts_sun_shade_sun/FigureS9_FittedModels_S.R"))

source(here("2_scripts_sun_shade_sun/FigureS10_FittedModels_Vcmaxt.R"))

source(here("2_scripts_sun_shade_sun/FigureS11_LightResponses.R"))

