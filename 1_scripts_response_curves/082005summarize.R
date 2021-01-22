#Processing cowpea gas exchange data for A/PPFD,  A/Ci,  and one point Vcmax
#light responses,  Aci responses and ci-correction of induction data
# all dealt with at the level of individual replicates
#082005nlme....R further process (or in the case of A/PPFD substitute for)
# these analyses prior to use in diurnal models (082005NaturePlantsdiurnal.R)

library(here)

#functions required for the analyses below

#AQ fitting
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/082005AQfit.R"))

#Acifit12.1
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/072020Acifit12.1.R"))

#Astar/one-point Vcmax
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/072020AstarVcmax.R"))

#specific to the analysis presented
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/functions_082005summarize.R"))

#Load and clean dataset (cp) and separate data objects for AC, AQ and inds
source(here("1_scripts_response_curves/082005cleanup.R"))

#a pdf that will hold all of the graphical output from analyses below
pdf(here("output/082005summarize.pdf"), h = 11, w = 8, paper = "a4")

#a list for outputs
cpGE <- list(NA)

################################################################################
################################################################################

#fit AQ responses

AQ.pars <- c("phi", "Asat", "theta", "Rd.AQ")

cpGE$AQ <- cpdf(AQ.pars)

######
#ALTHOUGH THE ABOVE IS MORE LONG-WIDED, IT REPLACES THE BELOW
#for (i in levels(cpGE$AQ$geno)){
#	par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
#
#	for (j in as.numeric(levels(factor(cpGE$AQ$block))) ){
#		AQin <- AQ[AQ$plant == paste(i, j, sep = "_"), ]
#		
#		if (nrow(AQin) > 1){
#			AQin <- AQin[c(1:14), ]#trim out the second 1500 when stomata are closed
#			m.AQin <- AQ.cost.fits(AQin,
#															params = c(phi = 0.05, Asat = 30, theta = 0.75, Rd = 1.9),
#															upp = NA,
#															low = NA,
#															Rd.fixed = F
#															)
#			cpGE$AQ[cpGE$AQ$plant == paste(i, j, sep = "_"), AQ.pars] <- m.AQin$coefs
#			plot.AQfit(m.AQin, xlim = c(0, 2000), ylim = c(-5, 45), phiPSII = TRUE)
#			title(main = paste(i, j, sep = "_"))
#			} else {
#					plot(1, 1, type = "n", axes = FALSE, bty = "n", main = j, xlab = "", ylab = "")
#					}
#}
#}

#cpGE$AQ

AQ.list <- by(AQ, AQ$plant, identity)

#first model
#for start values, see doAQfit in: functions_082005summarize.R
AQfits_no_fixed <- lapply(AQ.list, doAQfit)

#plot this
header.page(
  main = expression(italic(A)/PPFD~~responses~~all~~parameters~~estimated)
  )
plot.AQfit.bygeno(AQfits_no_fixed)

#summarize
cpGE$AQ <- mkAQsq(AQfits_no_fixed, AQ.pars)

################################################################################
################################################################################

#Establish mean values for operating states at start of A/ci

#calculate iWUE; why not?
AC$iWUE <- AC$A / AC$gsw

pvars <- c("A", "gsw", "Pci", "iWUE", "TleafCnd", "VPDleaf",
           "Fv.Fm", "PhiPS2", "NPQ", "qP_Fo", "qN_Fo"
)

cpGE$op <- cpdf(pvars)

cpGE$op[, pvars] <- data.frame(t(
  sapply(cpGE$op$plant, function(.){ AC[AC$plant == ., pvars][1, ] })
)
)

################################################################################
################################################################################
#fit Aci responses

#add a column to AC object that can be used to source Rd as wished
AC$Rd.AQ <- AC$plant
all(levels(AC$Rd.AQ) == cpGE$AQ$plant)
levels(AC$Rd.AQ) <- cpGE$AQ$Rd.AQ
AC$Rd.AQ
AC$Rd.AQ <- as.numeric(as.character(AC$Rd.AQ))
AC$Rd.AQ

#there are fewer plants than levels for plants
length(unique(AC$plant)) == length(levels(AC$plant))
unique(AC$plant) %in% levels(AC$plant)
levels(AC$plant) %in% unique(AC$plant)
#so, fix this
levels(AC$plant) <- replace(levels(AC$plant),
														!levels(AC$plant) %in% unique(AC$plant),
														NA
														)

#format AC as a list,  broken down by plant
AC.list <- by(AC, AC$plant, identity)

#first model
AC.fits_no_fixed <- lapply(AC.list, doACfit, 
														Rd.fixed = FALSE, 
														Gamma.star = NA,
														Gamma.star.fixed = FALSE, 
														Kco = NA,
														Kco.fixed = FALSE, 
														gm = NA,
														gm.fixed = FALSE
														)

#plot this
header.page(
  main = expression(
    italic(A)/italic(c)[i]~~responses~~all~~parameters~~estimated
    )
  )
plot.ACfit.bygeno(AC.fits_no_fixed)
dev.off()

#summarize
cpGE$AC_no_fixed <- mkACsq(AC.fits_no_fixed)

#next model
AC.fits_GammastarKco_fixed <- lapply(AC.list,
																			doACfit, 
																			Rd.fixed = FALSE, 
																			Gamma.star = NA,
																			Gamma.star.fixed = TRUE, 
																			Kco = NA,
																			Kco.fixed = TRUE, 
																			gm = NA,
																			gm.fixed = FALSE
																			)

header.page(
  main = expression(
  italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) *
    "* " * and~~italic(K)[CO]
  )
  )
plot.ACfit.bygeno(AC.fits_GammastarKco_fixed)

cpGE$AC_GammastarKco_fixed <- mkACsq(AC.fits_GammastarKco_fixed)

#next model
AC.fits_GammastarKcoRd_fixed <- lapply(AC.list, doACfit, 
																				Rd.fixed = TRUE, 
																				Gamma.star = NA,
																				Gamma.star.fixed = TRUE, 
																				Kco = NA,
																				Kco.fixed = TRUE, 
																				gm = NA,
																				gm.fixed = FALSE
																				)


header.page(
  main = expression(
  italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) *
    "*, " * ~~italic(K)[CO]~~and~~italic(R)[d]
  )
  )

plot.ACfit.bygeno(AC.fits_GammastarKcoRd_fixed)

cpGE$AC_GammastarKcoRd_fixed <- mkACsq(AC.fits_GammastarKcoRd_fixed)


#Asking: is gm functionally different between genotypes,
# or can it reasonably be set to an average value?

header.page(
  main = expression(
    atop("Is fixing " * italic(Gamma) * "* and" * italic(K)[CO] * " justified?",
         "Is it reasonable to fix " * italic(g)[m] * "?")
    )
  )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Kco ~ as.numeric(geno),
     data = cpGE$AC_no_fixed,
     log = "y",
     ylim = c(1, 110),
     las = 1,
     main = "no_fixed"
     )
plot(Kco ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     log = "y",
     ylim = c(1, 110),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Kco ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     log = "y",
     ylim = c(1, 110),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Gamma.star ~ as.numeric(geno),
     data = cpGE$AC_no_fixed,
     log = "y",
     ylim = c(0.1, 10),
     las = 1,
     main = "no_fixed"
     )
plot(Gamma.star ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     log = "y",
     ylim = c(0.1, 10),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Gamma.star ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     log = "y",
     ylim = c(0.1, 10),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_no_fixed,
     las = 1,
     main = "no_fixed"
     )
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )
plot(cpGE$AC_GammastarKco_fixed$Rd ~ cpGE$AC_GammastarKcoRd_fixed$Rd,
     ylim = c(0, 5),
     xlim = c(0, 5),
     las = 1
     )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_no_fixed,
     log = "y",
     ylim = c(5, 50),
     las = 1,
     main = "no_fixed"
     )
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     log = "y",
     ylim = c(5, 50),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     log = "y",
     ylim = c(5, 50),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(gm ~ as.numeric(geno),
     data = cpGE$AC_no_fixed,
     log = "y",
     las = 1,
     main = "no_fixed"
     )
plot(gm ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     log = "y",
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(gm ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     log = "y",
     las = 1,
     main = "GammastarKcoRd_fixed"
     )

#These plots indicate that Gammastar,  Kco,  and Rd
# are not identifiable by best fit from the data

#Fixing Gammastar and Kco enables Rd to be identified, 
# albeit with a slightly greater variance and greater top-end values.
#Notes:
#- fixing Rd sometimes results in poor fitting to the lowest ci data
#- When Kco and Gammastar are fixed, the estimate of Rd appears to fit data well

#It is probably reasonable to fix gm at a consistent value...

################################################################################
#Will fixing gm at a reasonable value improve the reliability of Rd,
# if Kco and Gammastar are also fixed?

summary(cpGE$AC_GammastarKcoRd_fixed$gm)
summary(cpGE$AC_GammastarKco_fixed$gm)
#a value of 5 mol/m2Pa is reasonable for gm,
# the median being 4.9-5 in constrained models
#(noting that the mean is strongly affected by outliers)

#mistkae identified here 0121 - Rd.fixed was set to TRUE, now corrected
AC.fits_GammastarKcogm_fixed <- lapply(AC.list,
																				doACfit, 
																				Rd.fixed = FALSE,																
																				Gamma.star = NA,
																				Gamma.star.fixed = TRUE, 
																				Kco = NA,
																				Kco.fixed = TRUE, 
																				gm = 5,
																				gm.fixed = TRUE
																				)

header.page(
  main = expression(
    italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) *
      " * , " * ~~italic(K)[CO]~~and~~italic(g)[m]
    )
  )

plot.ACfit.bygeno(AC.fits_GammastarKcogm_fixed)

cpGE$AC_GammastarKcogm_fixed <- mkACsq(AC.fits_GammastarKcogm_fixed)


#Asking whether fixing gm (Gammastar & Kco) results in reliable estimation of Rd
# as well as other params
header.page(main = expression("Impact of fixing " * italic(g)[m]))
	
par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
plot(cpGE$AC_GammastarKco_fixed$Rd ~ cpGE$AC_GammastarKcogm_fixed$Rd,
     ylim = c(0, 5),
     xlim = c(0, 5),
     las = 1
     )
#definitely does not introduce bias in Rd

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Gamma ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     ylim = c(0, 10),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Gamma ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     ylim = c(0, 10),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )
plot(Gamma ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 10),
     las = 1,
     main = "GammastarKcogm_fixed"
     )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     ylim = c(0, 50),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     ylim = c(0, 50),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 50),
     las = 1,
     main = "GammastarKcogm_fixed"
     )

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Vcmax ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(Vcmax ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcoRd_fixed"
     )
plot(Vcmax ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
#Vcmax estimates a bit more consistent with fixed gm

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(J ~ as.numeric(geno),
     data = cpGE$AC_GammastarKco_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKco_fixed"
     )
plot(J ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoRd_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcoRd_fixed")
plot(J ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
#no apparent impact on J

#As a sanity check,  since gm complicates the induction analysis,
# what happens if gm is fixed to ~Inf?
#0121 also spotted Rd = TRUE here, and corrected it 
AC.fits_GammastarKcoInfgm_fixed <- lapply(AC.list,
																					doACfit, 
																					Rd.fixed = FALSE, 
																					Gamma.star = NA,
																					Gamma.star.fixed = TRUE, 
																					Kco = NA,
																					Kco.fixed = TRUE, 
																					gm = 1e6,
																					gm.fixed = TRUE
																					)

header.page(
  main = expression(
    italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) *
      "*, " * ~~italic(K)[CO]~~and~~italic(infinite~~g)[m]
    )
  )

plot.ACfit.bygeno(AC.fits_GammastarKcoInfgm_fixed)

cpGE$AC_GammastarKcoInfgm_fixed <- mkACsq(AC.fits_GammastarKcoInfgm_fixed)

#right,  this has an impact!
table(cpGE$AC_GammastarKcogm_fixed$Ac.Aj >
        cpGE$AC_GammastarKcogm_fixed$Pci.op
      )
table(cpGE$AC_GammastarKcoInfgm_fixed$Ac.Aj >
        cpGE$AC_GammastarKcoInfgm_fixed$Pci.op
      )
#totally switches the limitation states...

header.page(main = expression("Impact of fixing " * italic(g)[m]~~to~~Inf))
	
par(mfrow = c(4, 2), mar = c(4, 5, 2, 1))
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
plot(Rd ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoInfgm_fixed,
     ylim = c(0, 5),
     las = 1,
     main = "GammastarKcoInfgm_fixed"
     )
points(Rd ~ as.numeric(geno),
       data = cpGE$AC_GammastarKcoRd_fixed,
       pch = 19,
       col = rgb(1, 0, 0, alpha = 0.3)
       )
legend(3.5, 4.5,
       xjust = 0.5,
       yjust = 0.5,
       bty = "n",
       pch = 19,
       col = rgb(1, 0, 0, alpha = 0.3),
       legend = expression(italic(R)[d]~~from~~light~~response)
       )

plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 50),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
plot(Ac.Aj ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoInfgm_fixed,
     ylim = c(0, 50),
     las = 1,
     main = "GammastarKcoInfgm_fixed"
     )

plot(Vcmax ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
plot(Vcmax ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoInfgm_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcoInfgm_fixed"
     )

plot(J ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcogm_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcogm_fixed"
     )
plot(J ~ as.numeric(geno),
     data = cpGE$AC_GammastarKcoInfgm_fixed,
     ylim = c(0, 400),
     las = 1,
     main = "GammastarKcoInfgm_fixed"
     )
#no apparent impact on J

#If gm = Inf,  can Kco, Gammastar & Rd all be estimated from the data?
AC.fits_Infgm_fixed <- lapply(AC.list,
															doACfit, 
															Rd.fixed = FALSE, 
															Gamma.star = NA,
															Gamma.star.fixed = FALSE, 
															Kco = NA,
															Kco.fixed = FALSE, 
															gm = 1e6,
															gm.fixed = TRUE
															)

header.page(
  main = expression(
    italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(infinite~~g)[m]
    )
  )

plot.ACfit.bygeno(AC.fits_Infgm_fixed)

cpGE$AC_Infgm_fixed <- mkACsq(AC.fits_GammastarKcoInfgm_fixed)

cpGE$AC_Infgm_fixed <- cpdf(AC.pars)

cpGE$AC_Infgm_fixed
#no,  Gammastar and Kco do behave better,  but Rd is collapsing

#Overall, for current purposes,
#cpGE$AC_GammastarKcoInfgm_fixed is a reliable description of the data.

#It implies Vcmax limitation in the steady state,
#so induction could be processed assuming minimal impacts of J

##############################
####I will need to choose which analysis I actually want
####Noting that there's currently a mistake in 082005NaturePlantsnlmeVcmax
#However, it would be wise to have an additional analysis including gm
# (noting that we don't really know how gm behaves during induction). 
#where since we don't know the actual dynamic behaviour of J and Vcmax but we can reasonably assume Vcmax recovers more slowly, 
#Vcmax half-time is estimated using only induction values where A<max(A[C]) from cpGE$AC_GammastarKcogm_fixed

################################################################################
################################################################################

#One-point Vcmax

#clean out Pci > 43 or < 0,  and A = 0
#hashed plots below were done in an interactive session
# removed here to keep them out of pdf
#plot(A ~ Pci, data = inds)
inds <- inds[inds$Pci < 44, ]
inds <- inds[inds$Pci > 0, ]
inds <- inds[inds$A > 0, ]
#plot(A ~ Pci, data = inds)
row.names(inds) <- c(1:nrow(inds))

#headers for calculated ouputs in induction dataset
Astar.pars <- c("nss.A", "nss.Pci", "Vcmax.t", "Astar")

#headers for new columns
#Vcmax.t  =  one-point estimate of Vcmax (based on re-arranged FvCB)
#Astar  =  A predicted at nss.Pci,  with Vcmax.t
#others are duplicates of info in data coded for easy plotting

#inputs
#list formatted induction data
inds[, Astar.pars] <- NA
inds.list <- by(inds, inds$plant, identity)

#list formatted Aci outputs
Aci_GammastarKcoInfgm.list <- by(cpGE$AC_GammastarKcoInfgm_fixed,
																	cpGE$AC_GammastarKcoInfgm_fixed$plant,
																	identity
																	)
#predict the values throughout induction for Inf gm model
inds.list_GammastarKcoInfgm_fixed <- mapply(do.Astarpredict,
																						inds.list,
																						Aci_GammastarKcoInfgm.list
																						)
#header for plot
header.page(
  main = expression(
    atop("induction with and without correction for " * italic(c)[i],
         with~~fixed~~italic(Gamma) * " * ,  " *
           ~~italic(K)[CO]~~and~~italic(infinite~~g)[m]
         )
    )
  )

#produce matching inputs
inds.list_GammastarKcoInfgm_fixed_paired <- trim.inds(
  inds.list_GammastarKcoInfgm_fixed, AC.fits_GammastarKcoInfgm_fixed
  )

AC.fits_GammastarKcoInfgm_fixed_paired <- trim.AC(
  AC.fits_GammastarKcoInfgm_fixed, inds.list_GammastarKcoInfgm_fixed
  )

#plot
plot.Astar.bygeno(inds.list_GammastarKcoInfgm_fixed_paired,
									AC.fits_GammastarKcoInfgm_fixed_paired
									)

#repeat for the fixed gm Aci
Aci_GammastarKcogm.list <- by(cpGE$AC_GammastarKcogm_fixed,
															cpGE$AC_GammastarKcogm_fixed$plant,
															identity
															)

inds.list_GammastarKcogm_fixed <- mapply(do.Astarpredict,
																					inds.list,
																					Aci_GammastarKcogm.list
																					)

header.page(
  main = expression(
    atop("induction with and without correction for " * italic(c)[i],
         with~~fixed~~italic(Gamma) * "*, " *
           ~~italic(K)[CO]~~and~~fixed~~finite~~italic(g)[m]
         )
    )
  )

inds.list_GammastarKcogm_fixed_paired <- trim.inds(
  inds.list_GammastarKcogm_fixed,	AC.fits_GammastarKcogm_fixed
  )

AC.fits_GammastarKcogm_fixed_paired <- trim.AC(
  AC.fits_GammastarKcogm_fixed, inds.list_GammastarKcogm_fixed
  )

plot.Astar.bygeno(inds.list_GammastarKcogm_fixed_paired,
									AC.fits_GammastarKcogm_fixed_paired
									)

#close the pdf
dev.off()

################################################################################
################################################################################

#the below go forwards for additional analysis in 082005NaturePlantsnlmeVcmax.R

save(inds.list_GammastarKcogm_fixed_paired,
			inds.list_GammastarKcoInfgm_fixed_paired,
			file = here("output/082005Astar.Rdata")
			)

save.image(here("output/082005summarize.Rdata"))
