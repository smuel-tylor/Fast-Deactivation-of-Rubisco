#Processing cowpea gas exchange data for A/PPFD,  A/Ci,  and one point Vcmax
#light responses,  Aci responses and ci-correction of induction data
# all dealt with at the level of individual replicates
#082005nlme....R further process (or in the case of A/PPFD substitute for)
# these analysesprior to use in diurnal models (082005NaturePlantsdiurnal.R)

library(here)

################################################################################
#Load and clean dataset (cp) and separate data objects for AC, AQ and inds
source(here("1_scripts_response_curves/082005cleanup.R"))

#a pdf that holds all of the graphical output from analyses below
pdf(here("output/082005summarize.pdf"), h = 11, w = 8, paper = "a4")

################################################################################
#establish mean values for operating states at start of A/ci

#calculate iWUE; why not?
AC$iWUE <- AC$A / AC$gsw

#generate a list that can be used to organise outputs
cpGE <- list(NA)

#within this,  a data frame for operating point values
names(cpGE) <- "op"

pvars <- c("A", "gsw", "Pci", "iWUE", "TleafCnd", "VPDleaf",
           "Fv.Fm", "PhiPS2", "NPQ", "qP_Fo", "qN_Fo"
           )

#function that generates dataframes for output
cpdf <- function(vars){
	fcs <- c("plant", "geno", "block", "cult")
	lp <- levels(cp$plant)
	tdf <- data.frame(
	  matrix(NA, length(lp), length(c(fcs, vars)))
	  )
	names(tdf) <- c(fcs, vars)
	tdf$plant <- lp
	geno.strings <- unlist(strsplit(lp, split = "_"))
	tdf$geno <- factor(geno.strings[seq(1, length(geno.strings), 2)],
	                   levels = levels(cp$geno)
	                   )
	tdf$block <- geno.strings[seq(2, length(geno.strings), 2)]
	tdf$cult <- replace(as.character(tdf$geno),
	                    tdf$geno == "Vadenantha"|tdf$geno == "TVNu-1948",
	                    "wild"
	                    )
	tdf$cult <- replace(tdf$cult, tdf$cult != "wild", "cultivar")
	tdf$cult <- factor(tdf$cult)
	tdf
}

cpGE$op <- cpdf(pvars)

cpGE$op[, pvars] <- data.frame(t(
  sapply(cpGE$op$plant, function(.){ AC[AC$plant == ., pvars][1, ] })
  )
  )

################################################################################
#fit AQ responses
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/082005AQfit.R"))

AQ.pars <- c("phi", "Asat", "theta", "Rd.AQ")

cpGE$AQ <- cpdf(AQ.pars)

doAQfit <- function(AQ, 
                    params = c(phi = 0.05, Asat = 30, theta = 0.75, Rd = 1.9),
                    upp = NA,
                    low = NA,
                    Rd.fixed = FALSE
){
  if (!is.null(AQ)){
    AQ.cost.fits(AQ,
                 params = params,
                 upp = upp,
                 low = low,
                 Rd.fixed = Rd.fixed
    )
  }
}

#pars here is the naming wanted for the coefs, if this differs from 
mkAQsq <- function(AQfits, pars = NA){
  
  #function to extract summary parameters from ACi fits into a dataframe row
  # and rename as necessary
  pull.AQpars <- function(AQfit, pars){
    if (!is.null(AQfit)){
      out <- AQfit$coefs
      if (all(!is.na(pars))){
        names(out) <- pars
      }
      out
    } else { rep(NA, 4) }
  }
  
  AQmat <- t(sapply(AQfits, pull.AQpars, pars = pars))
  AQdf <- data.frame(apply(AQmat, 2, unlist))
  AQdf$plant <- row.names(AQdf)
  merge(cpdf(NULL), AQdf, all.x = TRUE)
}

plotAQfit <- function(){
  
}

#function to plot Aci fits as one page per genotype,
# two panels per plant showing whole fit and intercept region
plot.AQfit2 <- function(AQfit){
  if (!is.null(AQfit)){
  plot.AQfit(AQfit,
             xlim = c(0, 1600),
             ylim = c(-5, 55),
             PhiPSII = TRUE
  )
  title(main = AQfit$data$plant[1])
  plot.AQfit(AQfit,
             xlim = c(0, 150),
             ylim = c(-5, 20),
             PhiPSII = FALSE
  )
  title(main = AQfit$data$plant[1])
  }
}

plot.AQfit.bygeno <- function(AQfits){
  for (i in levels(cp$geno)){
    par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
    lapply(AQfits[grep(i, names(AQfits), value = TRUE)], plot.AQfit2)
  }
}

#function to insert blank page in pdf output containing a title
# (main is an expression or character vector)
header.page <- function(main){
  par(mfrow = c(1, 1), mar = c(10, 10, 10, 10))
  plot(1, 1,
       type = "n",
       bty = "n",
       axes = FALSE,
       xlab = "",
       ylab = "",
       main = main
  )
}


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
AQfits_no_fixed <- lapply(AQ.list, doAQfit)

#plot this
header.page(
  main = expression(italic(A)/PPFD~~responses~~all~~parameters~~estimated)
  )
plot.AQfit.bygeno(AQfits_no_fixed)

#summarize
cpGE$AQ <- mkAQsq(AQfits_no_fixed, AQ.pars)

################################################################################
#fit Aci responses

#add a column to AC that can be used to source Rd as wished
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
#so
levels(AC$plant) <- replace(levels(AC$plant),
														!levels(AC$plant) %in% unique(AC$plant),
														NA
														)

#Code for Acifit12.1
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/072020Acifit12.1.R"))

#This function fits a curve using inputs and returns a fitted curve list object
#As above, the indexing and data selection are in a 'doACfits' function below
#List output means all elements of the fit are retained in these objects
#for subsequent plotting etc.
doACfit <- function(Aci, 
										Rd.fixed = FALSE, 
										Gamma.star = NA,
										Gamma.star.fixed = FALSE, 
										Kco = NA,
										Kco.fixed = FALSE, 
										gm = NA,
										gm.fixed = FALSE
										){

	A.op = Aci[1, "A"]
	Pci.op = Aci[1, "Pci"]
	Pca.op = Aci[1, "Pca"]
	
	Rd.guess <- Aci[1, "Rd.AQ"]
	if (is.na(Rd.guess)) {Rd.guess <- 2} #because Rd must be provided

	#generate the input
	#Gamma.star, Gamma.star.fixed, Kco, Kco.fixed, gm, gm.fixed, var.J, stoich
	input.AC <- Acifit.input(Aci, 
														Rd = Rd.guess,
														Rd.fixed, 
														Gamma.star,
														Gamma.star.fixed, 
														Kco,
														Kco.fixed, 
														gm,
														gm.fixed, 
														stoich = "NADPH", 
														var.J = FALSE
														)
	#add in lim states
	input.AC <- gen.lim.states(input.AC)
	#cost fits
	input.AC <- Aci.cost.fits(input.AC)
	#admissibility
	input.AC <- Acifit.admit(input.AC)
	#colimited points
	input.AC <- Acifit.colimit(input.AC)
	#choose model
	AdT <- input.AC$models[input.AC$models$admit == TRUE, ]
	best.AC = AdT[AdT$cost == min(AdT$cost, na.rm = TRUE), ]
	input.AC <- Acifit.evaluate(best.AC, input.AC, A.op, Pci.op, Pca.op)
	
	input.AC
}

#function to make a summary of key params from fits... complicated,
# but necessary to make a square df
mkACsq <- function(ACfits){
	
	#parameters for Aci fits that are useful for summaries
	AC.pars = c("Vcmax", "J", "TPU", "Rd", "Gamma.star", "Kco", "gm",
	            "A.40", "A.op", "Ls", "Pci.op",
	            "Gamma", "Ac.Aj", "Aj.Ap", "Ac.Ap",
	            "Vcmax.25", "J.25", "TPU.25", "Rd.25", "gm.25"
	            )
	
	#extract summary parameters from ACi fits into a dataframe row
	pull.ACpars <- function(input.AC, pars){
	  core.pars <- c("Vcmax", "J", "TPU", "Rd", "Gamma.star", "Kco", "gm")
		out <- c(input.AC$chosen.mod[ , core.pars], 
							input.AC$F_Slimitation, 
							input.AC$ci.transitions, 
							input.AC$temperature.normalised
							)
		names(out) <- pars
		out
	}
	
	ACmat <- t(sapply(ACfits, pull.ACpars, pars = AC.pars))
	ACdf <- data.frame(apply(ACmat, 2, unlist))
	ACdf$plant <- row.names(ACdf)
	merge(cpdf(NULL), ACdf, all.x = TRUE)
}

#function to plot Aci fits as one page per genotype,  two panels per plant
# showing whole fit and Gamma region
plot.ACfit2 <- function(input.AC){

	plot.Acifit(input.AC,
	            xlim = c(0, 90),
	            ylim = c(-5, 55),
	            PhiPSII = TRUE,
	            F_Slimitation = TRUE
	            )
  title(main = input.AC$data$plant[1])
	plot.Acifit(input.AC,
	            xlim = c(0, 20),
	            ylim = c(-5, 20),
	            PhiPSII = FALSE,
	            F_Slimitation = FALSE
	            )
	title(main = input.AC$data$plant[1])

	}

plot.ACfit.bygeno <- function(AC.fits){
	for (i in levels(cp$geno)){
		par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
		lapply(AC.fits[grep(i, names(AC.fits), value = TRUE)], plot.ACfit2)
	}
}
	
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
  main = expression(italic(A)/italic(c)[i]~~responses~~all~~parameters~~estimated)
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

header.page(main = expression(italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) * "* " * and~~italic(K)[CO]))
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


header.page(main = expression(italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) * "*, " * ~~italic(K)[CO]~~and~~italic(R)[d]))
plot.ACfit.bygeno(AC.fits_GammastarKcoRd_fixed)

cpGE$AC_GammastarKcoRd_fixed <- mkACsq(AC.fits_GammastarKcoRd_fixed)


#is gm functionally different between genotypes,  or can it reasonably be set to an average value?

header.page(main = expression(atop("Is fixing " * italic(Gamma) * "* and" * italic(K)[CO] * " justified?", "Is it reasonable to fix " * italic(g)[m] * "?")))
	
par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Kco~as.numeric(geno), data = cpGE$AC_no_fixed, log = "y", ylim = c(1, 110), las = 1, main = "no_fixed")
plot(Kco~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, log = "y", ylim = c(1, 110), las = 1, main = "GammastarKco_fixed")
plot(Kco~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, log = "y", ylim = c(1, 110), las = 1, main = "GammastarKcoRd_fixed")

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Gamma.star~as.numeric(geno), data = cpGE$AC_no_fixed, log = "y", ylim = c(0.1, 10), las = 1, main = "no_fixed")
plot(Gamma.star~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, log = "y", ylim = c(0.1, 10), las = 1, main = "GammastarKco_fixed")
plot(Gamma.star~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, log = "y", ylim = c(0.1, 10), las = 1, main = "GammastarKcoRd_fixed")

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Rd~as.numeric(geno), data = cpGE$AC_no_fixed, las = 1, main = "no_fixed")
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, ylim = c(0, 5), las = 1, main = "GammastarKco_fixed")
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, ylim = c(0, 5), las = 1, main = "GammastarKcoRd_fixed")
plot(cpGE$AC_GammastarKco_fixed$Rd ~ cpGE$AC_GammastarKcoRd_fixed$Rd, ylim = c(0, 5), xlim = c(0, 5), las = 1)

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_no_fixed, log = "y", ylim = c(5, 50), las = 1, main = "no_fixed")
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, log = "y", ylim = c(5, 50), las = 1, main = "GammastarKco_fixed")
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, log = "y", ylim = c(5, 50), las = 1, main = "GammastarKcoRd_fixed")

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(gm~as.numeric(geno), data = cpGE$AC_no_fixed, log = "y", las = 1, main = "no_fixed")
plot(gm~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, log = "y", las = 1, main = "GammastarKco_fixed")
plot(gm~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, log = "y", las = 1, main = "GammastarKcoRd_fixed")

#these plots indicate that Gammastar,  Kco,  and Rd are not identifiable by best fit from the data

#Fixing Gammastar and Kco enables Rd to be identified, 
#albeit with a slightly greater variance and greater top-end values
#Note that Fixing Rd sometimes results in poor fitting to the lowest point on the A/ci response
#the higher estimate of Rd when Kco and Gammastar are fixed is a good representation of the data

#It is probably reasonable to fix gm at a consistent value...
#Will this improve the reliability of,  e.g.,  Rd if Kco and Gammastar are also fixed?

summary(cpGE$AC_GammastarKcoRd_fixed$gm)
summary(cpGE$AC_GammastarKco_fixed$gm)
#a value of 5 mol/m2Pa is reasonable for gm,  as the median is 4.9-5 in constrained models (note the mean is strongly affected by rare outliers

AC.fits_GammastarKcogm_fixed <- lapply(AC.list,
																				doACfit, 
																				Rd.fixed = TRUE,																
																				Gamma.star = NA,
																				Gamma.star.fixed = TRUE, 
																				Kco = NA,
																				Kco.fixed = TRUE, 
																				gm = 5,
																				gm.fixed = TRUE
																				)

header.page(main = expression(italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) * " * , " * ~~italic(K)[CO]~~and~~italic(g)[m]))
plot.ACfit.bygeno(AC.fits_GammastarKcogm_fixed)

cpGE$AC_GammastarKcogm_fixed <- mkACsq(AC.fits_GammastarKcogm_fixed)


#####
header.page(main = expression("Impact of fixing " * italic(g)[m]))
	
par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, ylim = c(0, 5), las = 1, main = "GammastarKco_fixed")
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, ylim = c(0, 5), las = 1, main = "GammastarKcoRd_fixed")
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 5), las = 1, main = "GammastarKcogm_fixed")
plot(cpGE$AC_GammastarKco_fixed$Rd ~ cpGE$AC_GammastarKcogm_fixed$Rd, ylim = c(0, 5), xlim = c(0, 5), las = 1)
#definitely does not introduce bias in Rd

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Gamma~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, ylim = c(0, 10), las = 1, main = "GammastarKco_fixed")
plot(Gamma~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, ylim = c(0, 10), las = 1, main = "GammastarKcoRd_fixed")
plot(Gamma~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 10), las = 1, main = "GammastarKcogm_fixed")

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, ylim = c(0, 50), las = 1, main = "GammastarKco_fixed")
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, ylim = c(0, 50), las = 1, main = "GammastarKcoRd_fixed")
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 50), las = 1, main = "GammastarKcogm_fixed")

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(Vcmax~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, ylim = c(0, 400), las = 1, main = "GammastarKco_fixed")
plot(Vcmax~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcoRd_fixed")
plot(Vcmax~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcogm_fixed")
#Vcmax estimates a bit more consistent with fixed gm

par(mfrow = c(3, 2), mar = c(4, 5, 2, 1))
plot(J~as.numeric(geno), data = cpGE$AC_GammastarKco_fixed, ylim = c(0, 400), las = 1, main = "GammastarKco_fixed")
plot(J~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcoRd_fixed")
plot(J~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcogm_fixed")
#no apparent impact on J

#As a sanity check,  since gm complicates the induction analysis,  what happens if gm is fixed to ~Inf
AC.fits_GammastarKcoInfgm_fixed <- lapply(AC.list,
																					doACfit, 
																					Rd.fixed = TRUE, 
																					Gamma.star = NA,
																					Gamma.star.fixed = TRUE, 
																					Kco = NA,
																					Kco.fixed = TRUE, 
																					gm = 1e6,
																					gm.fixed = TRUE
																					)

header.page(main = expression(italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(Gamma) * "*, " * ~~italic(K)[CO]~~and~~italic(infinite~~g)[m]))
plot.ACfit.bygeno(AC.fits_GammastarKcoInfgm_fixed)

cpGE$AC_GammastarKcoInfgm_fixed <- mkACsq(AC.fits_GammastarKcoInfgm_fixed)

#right,  this has an impact!
table(cpGE$AC_GammastarKcogm_fixed$Ac.Aj > cpGE$AC_GammastarKcogm_fixed$Pci.op)
table(cpGE$AC_GammastarKcoInfgm_fixed$Ac.Aj > cpGE$AC_GammastarKcoInfgm_fixed$Pci.op)
#totally switches the limitation states...

header.page(main = expression("Impact of fixing " * italic(g)[m]~~to~~Inf))
	
par(mfrow = c(4, 2), mar = c(4, 5, 2, 1))
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 5), las = 1, main = "GammastarKcogm_fixed")
plot(Rd~as.numeric(geno), data = cpGE$AC_GammastarKcoInfgm_fixed, ylim = c(0, 5), las = 1, main = "GammastarKcoInfgm_fixed")
points(Rd~as.numeric(geno), data = cpGE$AC_GammastarKcoRd_fixed, pch = 19, col = rgb(1, 0, 0, alpha = 0.3))
legend(3.5, 4.5, xjust = 0.5, yjust = 0.5, bty = "n", pch = 19, col = rgb(1, 0, 0, alpha = 0.3), legend = expression(italic(R)[d]~~from~~light~~response))

plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 50), las = 1, main = "GammastarKcogm_fixed")
plot(Ac.Aj~as.numeric(geno), data = cpGE$AC_GammastarKcoInfgm_fixed, ylim = c(0, 50), las = 1, main = "GammastarKcoInfgm_fixed")

plot(Vcmax~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcogm_fixed")
plot(Vcmax~as.numeric(geno), data = cpGE$AC_GammastarKcoInfgm_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcoInfgm_fixed")

plot(J~as.numeric(geno), data = cpGE$AC_GammastarKcogm_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcogm_fixed")
plot(J~as.numeric(geno), data = cpGE$AC_GammastarKcoInfgm_fixed, ylim = c(0, 400), las = 1, main = "GammastarKcoInfgm_fixed")
#no apparent impact on J

#If gm is assumed infinite,  can Kco,  Gammastar & Rd all be estimated from the data?
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

header.page(main = expression(italic(A)/italic(c)[i]~~responses~~with~~fixed~~italic(infinite~~g)[m]))
plot.ACfit.bygeno(AC.fits_Infgm_fixed)

cpGE$AC_Infgm_fixed <- mkACsq(AC.fits_GammastarKcoInfgm_fixed)

cpGE$AC_Infgm_fixed <- cpdf(AC.pars)

cpGE$AC_Infgm_fixed
#no,  Gammastar and Kco do behave better,  but Rd is collapsing

#overall,  cpGE$AC_GammastarKcoInfgm_fixed is a reliable description of the data
#it implies Vcmax limitation in the steady state,  so induction could be processed assuming minimal impacts of J
#However,  it would be wise to have an additional analysis including gm (noting that we don't really know how gm behaves in these circumstances), 
#where since we don't know the actual dynamic behaviour of J and Vcmax but we can reasonably assume Vcmax recovers more slowly, 
#Vcmax half-time is estimated using only induction values where A<max(A[C]) from cpGE$AC_GammastarKcogm_fixed

#develop Astar for these two analyses using one-point Vcmax
source("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rscripts/072020AstarVcmax.R")

#dataframe for induction summaries
inds <- cp[cp$curve == "induction", ]
#clean out Pci>43 or <0,  and A = 0
#rmoving plots below because no longer relying on these and don't need them in the pdf output
#plot(A~Pci, data = inds)
inds <- inds[inds$Pci < 44, ]
inds <- inds[inds$Pci > 0, ]
inds <- inds[inds$A > 0, ]
#plot(A~Pci, data = inds)
row.names(inds) <- c(1:nrow(inds))

#headers for calculated ouputs in induction dataset
Astar.pars <- c("nss.A", "nss.Pci", "Vcmax.t", "Astar")
#headers for new columns
#Vcmax.t  =  one-point estimate of Vcmax (based on re-arranged FvCB)
#Astar  =  A predicted at nss.Pci,  with Vcmax.t
#others are duplicates of info in data coded for easy plotting

#function that adds Astar and associated values to dataframe
do.Astarpredict <- function(inds, Aci){
	
	Aci.mod = Aci
	stoich = "NADPH"
	A = inds$A
	Pci = inds$Pci
	Pci.ss = Aci$Pci.op
	
	inds[, Astar.pars] <- one.point.Vcmax(Aci.mod, stoich, A, Pci, Pci.ss)
	inds
}

#plotting functions
plot.Astar2 <- function(inds, input.AC){

	plot.Acifit(input.AC,
							xlim = c(0, 90),
							ylim = c(-5, 55),
							phiPSII = TRUE,
							F_Slimitation = TRUE
							)
	points(A~Pci, data = inds, col = 3, pch = 21)
	title(main = input.AC$data$plant[1])
	
	plot(A~induction.s,
				data = inds, 
				xlim = c(-60, 900),
				ylim = c(0, 40), 
				xlab = expression(time~~since~~shade~~(s)), 
				ylab = expression(italic(A)~~(mu * mol~~m^-2~~s^-1)), 
				las = 1
				)
	points(Astar~induction.s, data = inds, col = 2, pch = 21)
	legend(900, 0, bty = "n", xjust = 1, yjust = 0, pch = rep(21, 2), col = c(1, 2),
					legend = expression(italic(A), italic(c)[i] * " corrected,  " * italic(A) * " * ")
					)
	title(main = input.AC$data$plant[1])

	}

plot.Astar.bygeno <- function(inds.list, input.AC.list){
	for (i in levels(cp$geno)){
		par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
		mapply(plot.Astar2,
						inds = inds.list[grep(i, names(inds.list), value = TRUE)],
						input.AC = input.AC.list[grep(i, names(input.AC.list), value = TRUE)]
						)
	}
}

#to make the above plot function work,  find the sets that have both an induction and an Aci response fitted
trim.inds <- function(inds.list, input.AC.list){
	if ( length(input.AC.list) < length(inds.list) ){
		inds.list <- inds.list[names(input.AC.list)]
	}
		
	if ( length(input.AC.list) > length(inds.list) ){
		input.AC.list <- input.AC.list[names(inds.list)]
	}
	
	inll <- !unlist(lapply(inds.list, is.null))
	anll <- !unlist(lapply(input.AC.list, is.null))
	
	inds.list[inll & anll]

}

trim.AC <- function(input.AC.list, inds.list){
	if ( length(inds.list) < length(input.AC.list) ){
		input.AC.list <- input.AC.list[names(inds.AC.list)]
	}
		
	if ( length(inds.list) > length(input.AC.list) ){
		inds.list <- inds.list[names(input.AC.list)]
	}
	
	inll <- !unlist(lapply(inds.list, is.null))
	anll <- !unlist(lapply(input.AC.list, is.null))
	
	input.AC.list[inll & anll]
}

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
header.page(main = expression(
															atop("induction with and without correction for " * italic(c)[i], 
																		with~~fixed~~italic(Gamma) * " * ,  " * ~~italic(K)[CO]~~and~~italic(infinite~~g)[m]
																		)
															)
						)
#produce matching inputs
inds.list_GammastarKcoInfgm_fixed_paired <- trim.inds(inds.list_GammastarKcoInfgm_fixed,
																											AC.fits_GammastarKcoInfgm_fixed
																											)
AC.fits_GammastarKcoInfgm_fixed_paired <- trim.AC(AC.fits_GammastarKcoInfgm_fixed,
																									inds.list_GammastarKcoInfgm_fixed
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

header.page(main = expression(
															atop("induction with and without correction for " * italic(c)[i],
																		with~~fixed~~italic(Gamma) * "*, " * ~~italic(K)[CO]~~and~~fixed~~finite~~italic(g)[m]
																		)
															)
						)

inds.list_GammastarKcogm_fixed_paired <- trim.inds(inds.list_GammastarKcogm_fixed,
																										AC.fits_GammastarKcogm_fixed
																										)
AC.fits_GammastarKcogm_fixed_paired <- trim.AC(AC.fits_GammastarKcogm_fixed,
																								inds.list_GammastarKcogm_fixed
																								)

plot.Astar.bygeno(inds.list_GammastarKcogm_fixed_paired,
									AC.fits_GammastarKcogm_fixed_paired
									)


dev.off()
#these are taken forwards for additional analysis in 072006timesequencemodel_gasexchange.R

save(inds.list_GammastarKcogm_fixed_paired,
			inds.list_GammastarKcoInfgm_fixed_paired,
			file = "082005Astar.Rdata"
			)

save.image("082005summarize.Rdata")
