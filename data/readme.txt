Data used in the Taylor et al Nature Plants Rubisco Activation analysis

Data_FigS2.csv
See also the published paper and supplementary materials
Sample - indicates one of three treatments:
  - MES - leaf discs cut from a leaf and incubated under greenhouse lights
          on the surface of a beaker of MES buffer for 60 min before freezing
  - WATER - leaf discs cut from a leaf and incubated under greenhouse lights
            on the surface of a beaker of water for 60 min before freezing 
  - FROZEN - leaf discs cut from a leaf and snap frozen directly  
Genotype - a single genotype IT86D-1010 was used in this experiment
Vi - initial activity of Rubisco (micromol/m2s)
Vt - total activity of Rubisco (micromol/m2s)
Act - Vi/Vt 
Act.S - Act * 100

Data_FigS3.csv
Condition - duration over which leaf discs were incubated on a MES buffer under
            lights
ActS - Vi_umol_m2_s/Vt3_umol_m2_s
Activation State - ActS * 100
Vi_umol_min_mL - initial activity of Rubisco expressed per extraction volume
Vt_umol_min_mL - total activity of Rubisco expressed per extraction volume
Vi_umol_min_mg - initial activity of Rubisco per total soluble protein
Vt_umol_min_mg - total activity of Rubisco per total soluble protein
Vi_umol_m2_s - initial activity of Rubisco expressed on a leaf area basis
Vt3_umol_m2_s - total activity of Rubisco expressed on a leaf area basis
                - 3 indicates duration of incubation before quenching of reaction

Data_FigS4.csv
Genotype - four genotypes were used:
            - IT86D-1010
            - IT82E-16
            - TVNu-1948 (aka V. Sp Savi)
            - V. adenantha
Plant - code indicating which plant within each genotype was sampled,
        synonymous with leaf
Replicate - four replicate discs were cut per plant leaf and treatment
Bicarbonate_mM - the treatment of interest: concentration of bicarbonate added
                  to the incubation buffer
Duration_min - length of incubation
Vi_umol_min_mL - initial activity based on extract volume
Vt_umol_min_mL - total activity based on extract volume
ActivationState - initial/total activity * 100

Data_FigS5_IT82E16.csv
Data_FigS5_IT86D1010.csv
Data_FigS5_TVNu1948.csv
Data_FigS5_Vadenantha.csv
Extract - ID attributed the the individual extraction
SampleID - ID attributed to the sample from which the extraction was prepared
Genotype - aligns with filename, one genotype per file
Plant - code indicating which plant in the greenhouse was sampled,
        synonymous with leaf 
Replicate - code identifying individual discs cut from a single Plant
Treatment - light level, measured at the location where the leaf disc was
            subsequently incubated while floating on a buffer (micromol/m2s) 
Vi_umol_min_mL - initial activity based on extract volume
Vt_umol_min_mL - total activity based on extract volume
Vi_umol_min_mg - initial activity per mass of total soluble protein
Vt3_umol_min_mg - total activity per mass of total soluble protein
Vi_umol_m2_s - initial activity of Rubisco expressed on a leaf area basis
Vt3_umol_m2_s - total activity of Rubisco expressed on a leaf area basis
ActS_ratio - Vi/Vt
ActS - ActS_ratio * 100

level2.csv
Ray tracing model output for second layer of soybean canopy from
  Zhu et al 2004 J Exp Bot 55:1167
Time.h - time of day in hours
PPFD - micromol/m2s

RIPE_20210118_cowpeaGE.csv
Compiled gas exchange data used in analyses. Most common columns from LI-6800F
  output are included but system constants from file headers are not included.
  For many column headers, the best reference may be system documentation.
  LI-6800F software was:
  - Bluestem v.1.2.2
  - Scripts version 2017.12 1.2.1, Oct 2017
  - Fluorometer version was 1.1.6
comments
geno - cowpea accession, one of:- IT86D-1010
                                - TVNu-1948
                                - V. adenantha
                                - IT82E-16
curve - describes which response curve was being measured:
          - induction - sun-shade-sun assay
          - AQ - PPFD response
          - Aci - CO2 response
block - identifies a batch of plants grown in the glasshouse
obs - LI-6800F file observation number
time - time from LI-6800F
elapsed - elapsed time since LI-6800F file was opened
induction.s - for curve = induction, timepoint in s from end of shade during
              the sun-shade sun sequence
Qin - PPFD in cuvette micromol/m2s
date - from LI-6800F OS
TIME - from LI-6800F OS
E	- transpiration mmol/m2s
A	- net CO2 assimilation micromol/m2s
Ca	- cuvette [CO2] micromol/mol
Ci	- leaf intercellular [CO2] micromol/mol
Pci	- leaf intercellular [CO2] Pa
Pca	- cuvette [CO2] Pa
gsw - stomatal conductance to H2O (mol/m2s)
gbw - boundary layer conductance to H2O
gtw - total conductance to H2O
gtc - total conductance to CO2
Rabs - energy absorption by leaf W/m2
TleafEB - Leaf temperature from energy balance, here equiv to TleafCnd (Celcius)
TleafCnd - LEaf temperature by thermocouple contact (Celcius)
SVPleaf - saturated vapour pressure for the leaf interior (kPa)
RHcham - relative humidity in the leaf cuvette (%)
VPcham - vapour pressure in the leaf cuvette (kPa)
SVPcham -  saturated vapour pressure for the leaf interior (kPa)
VPDleaf - vapour pressure deficit from leaf interior to air (kPa)
Leak - LI-6800F variable indicating amount of outwards leakeage form leaf cuvette
LeakPct - LI-6800F variable for Leak, expressed as percent
CorrFact - LI-6800F variable used for leak correction
CorrFactPct - LI-6800F variable used for leak correction, expressed as percent
Fan - LI-6800F mixing fan output (micromol/s)
DarkAdaptedID - Identifier for Fo, Fm, Fv measure being used in calculations
Fo - dark adapted fluorescence
Fm - maximum fluorescence during a saturating flash, dark adapted
Fv - Fm - Fo
Fv/Fm - as stated
LightAdaptedID - Identifier for Fs and Fm' measures being used in calculations
Fs - fluorescence in the light
Fm' - maximum fluorescence during a saturating flash, light adapted
PhiPS2 - light adapted quantum yield of PSII
PS2/1 - used in calculations of ETR
Qabs_fs - not used here, PPFD absorbed at Fs (micromol/m2s)
Afs - Net CO2 assimilation at Fs (micromol/m2s)
ETR - electron transport rate (micromol/m2s)
Fv'/Fm' - (Fm' - Fo') / Fm' (note Fo' not measured here)
NPQ - (Fm - Fm') / Fm
qP_Fo - (Fm' - Fs) / (Fm' - Fo)
qN_Fo - (Fm - Fm') / (Fm - Fo)
Qabs - LI-6800F approximation of leaf absorbed PPFD (micromol/m2s)
alpha - coefficient used to calculate Qabs
convert - used to convert between PPFD and energy (J/micromol)
TIME - LI-6800F system variable
CO2_s - [CO2] sample IRGA (micromol/mol)
CO2_r - [CO2] reference IRGA (micromol/mol)
H2O_s - [H2O] sample IRGA (mmol/mol)
H2O_r	- [H2O] reference IRGA (mmol/mol)
Flow - rate of airflow through system (micromol/s)
Pa - Air Pressure measured in system (kPa)
dPcham	- positive pressude differential in cuvette (kPa)
Tair - Temperature of air in IRGAs (Celcius)
Tleaf - Temperature of leaf, synonymous with TleafCnd (Celcius)
Fan_speed	- chamber fan speed (rpm)
Qamb_out - PPFD at sensor mounted on exterior of head (micromol/m2s)
Î”CO2 - DeltaCO2, differential between reference and sample (micromol/mol)
CO2_s_d - sample [CO2] used for DeltaCO2 (micromol/mol)
CO2_r_d - reference [CO2] used for DeltaCO2 (micromol/mol)
Î”H2O - DeltaH2O, differential between reference and sample
[CO2_a, CO2_b, H2O_a, H2O_b, e_s, e_r, Td_s, Td_r, Q, f_red, f_blue, f_farred]
- LI-6800F system
F - fluorescence
Q_modavg - average PPFD from modulated light used in fluorescence (micromol/m2s)
[F_dc, Pc, Tled, TDigital, TPreamp, TPwrSpy, TDrive, Q_red, Q_blue, Q_farred,
TSPF, F_avg, dF/dt, dF_dc/dt, F_dc_avg, period, time, hhmmss, count, co2_adj,
h2o_adj, co2_at, h2o_at, co2_cv, h2o_cv, CO2_r:MN, CO2_r:SLP, CO2_r:SD,
CO2_r:OK,	CO2_s:MN, CO2_s:SLP, CO2_s:SD, CO2_s:OK, H2O_r:MN, H2O_r:SLP,
H2O_r:SD, H2O_r:OK,	H2O_s:MN, H2O_s:SLP, H2O_s:SD, H2O_s:OK, Stable, Total,
State]
- LI-6800F system & stability
MatchValveR - match valve status Reference (%)
MatchValveS - match valve status Sample (%)
MatchCO2 - CO2 match offset (micomol/mol)
MatchH2O - H2O match offset (mmol/mol)
DIAG - match valve position
Flow_s - sample IRGA flow (micromol/s)
Flow_r - reference IRGA flow  (micromol/s)
Txchg - tmeperature exchanger temperature (Celcius)
Tirga - gas analyser temperature (Celcius)
[Tchopper,	Ts, Tr, Txchg_sp, CO2_r_sp, H2O_r_sp, SS_s , SS_r]
- LI-6800F system

RIPE_20210121_cowpeaAS.csv
Complete Rubisco Activation time series dataset for the sun-shade-sun experiment
  described in the manuscript
comments - a note on file history
extract - code given to extraction
sampleID	- Identifier
geno - cowpea accession, one of:- IT86D-1010
                                - TVNu-1948
                                - V. adenantha
                                - IT82E-16
block - identifies a plant within the greenhouse
plant - idnetifier combining geno and block
time - point during the sun-shade-sun sequence (min)
Vi_umol_min_mL - initial activity based on extract volume
Vt3_umol_min_mL - initial activity based on extract volume
ActS_ratio - Vi_umol_min_mL/Vt3_umol_min_mL
ActivationState - ActS_ratio * 100 
