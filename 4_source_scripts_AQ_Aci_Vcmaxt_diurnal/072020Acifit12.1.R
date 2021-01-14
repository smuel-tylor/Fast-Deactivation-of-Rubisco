#Functions for Aci fitting with Gu optimisation approach (Gu et al., 2010 PC&E)
#using TPU assuming alpha=0

################################################################################ 
##notes from previous version (081913Acifit12.1.R)
##
##- Redeveloped from previous versions by breaking down further into functions
##
##- Allows A.op to be optional
##  for use with dynamic A/ci data for which A.op is not appropriate
##
##- Test for admissibility uses carboxylation states
##  and fixed to prevent issues with admissibility == NA if all params estimated
##
##- Allows flexibility in parameter choices
##  equations modified for ATP rather than NADPH stoichiometry
##  you can fit or estimate:
##  Rd
##  Gamma.star
##  Kco
##  (choice of these impacts N per segment)
##  gm - if specified as NULL = fit Rd and use Bernacchi gm after Sharkey 2007
##
##- Uses partial pressures
##
##- Allows any combination of limitation states, including lack of one or two
##
##- Has optional plotting
##  Farquhar Sharkey limitation plotting fixed 0319
##  updated June 2019 to allow legend to be an option and include ...
##  modified 0819 to enable co-limited points to be fit
##  (this is annotated: use find 0819 to identify modified segments)
################################################################################

################################################################################
##Other, more recent changes
##032027 Fixed Acifit.evaluate
##  set A.op and Pci.op defaults to NULL, outputting NA for Ls in these cases
##
##072020
##- Further updates on Acifit.evaluate and plot.Acifit
##  created A.ca, replacing what had been A.40
##  resulting in more correct Farquahar-Sharkey fits
##
##- Added coef.set option to Acifit.input
##  lets you choose coefficients: default "Sharkey", or "Silva-Perez" (wheat)
##
##- Fixed an issue with Acifit.input Pair
##  in earlier version this was always derived from Pci and Ci
##
##- 0121 Updated style (though not code itself)
## for easier reading, after submission to Nature Plants
################################################################################


#A function that sets up the inputs needed for Aci analysis
#doing basic unit checking and identifying options

#requires a data frame with
# A (micromol m^-2 s^-1)
# Ci (Pa) as Pci
# Pair (kPa)
# Tleaf (degrees C)
#If variable J gm is to be estimated (not tested recently/in this version)
# ETR (micromol m^-2 s^-1)
#For validation, phiPSII may also be supplied
# this can be used in plotting functions but is not used analytically
#Providing Gamma.star or Kco as NULL will fix these as per coef.set
#Providing any values as NA will result in parameters being estimated
# using coef.set estimate as a starting point
#Providing values will result in parameter being estimated
# using the provided value as a starting point, unless parameter.fixed = TRUE
#Providing Rd or gm as NULL will result in these being estimated
# using 'reasonable' starting values
#If var.J = TRUE
# an attempt is made to estimate gm using variable J (not tested/validated)


Acifit.input <- function(data,
                         Rd,
                         Rd.fixed,
                         Gamma.star,
                         Gamma.star.fixed,
                         Kco,
                         Kco.fixed,
                         gm,
                         gm.fixed,
                         var.J,
                         stoich,
                         coef.set = "Sharkey"
                         ){
		
	#Stoichiometry
  #defaults to ATP because it was used in the original Bernacchi derivation
  #of temperature functions for J
	if (any(is.null(stoich), is.na(stoich))) { stoich = "ATP" } else {
	                                                              stoich = stoich}
	
	#basic parameter set that will always require starting values
	params = list(Vcmax = 200, J = 250, TPU = 12)
	
	#if any values are to be fixed, set up a list of fixed values
	if (any(Rd.fixed == TRUE,
	        Gamma.star.fixed == TRUE,
	        is.null(Gamma.star),
	        Kco.fixed == TRUE,
	        is.null(Kco),
	        gm.fixed == T
	        )) { fixed = list(NA) } else { fixed = NULL }
		
	#format data
	Acidat <- data #should include:A, Pci, Tleaf, Pair as columns
	Acidat <- Acidat[order(Acidat[ , "Pci"]), ] #ensure order by Ci

	#if Pair is not in the file it is Pci / (Ci * 1e-6) Pa
	#correction made here in 0720 rewrite from any -> all
	if (all(names(Acidat) != "Pair")){
	Acidat$Pair <- rep(NA, nrow(Acidat))
	Acidat$Pair <- Acidat$Pci / (Acidat$Ci * 1e-6) * 1e-3 #kPa match LI-COR output
	}
	#print(Acidat)
	#plot(A~Ci,data=Acidat)
	
	##physical constants for fit
	P = mean(Acidat[ , "Pair"]) * 1e3 #kPa to Pa
	#note that in new LICOR output Pair would be the sum of Pair and deltaP
	#the latter is usually 0.1-0.2 kPa
	TleafK = mean(Acidat[ , "Tleaf"]) + 273.15 #mean leaf temperature in Kelvin
	O2 = 209.5 * 1e-3 * P #mmol / mol -> mol / mol -> Pa
	R = 8.314E-3 #gas constant kJ / molK
	physical = list(P.Pa = P, Tleaf.K = TleafK, O2 = O2, R = R)
	#print(physical)
	
	##fixed parameters after Bernacchi 2001, from Sharkey et al 2007 PC&E 
	##fitting to pa CO2 Cc
	
	#Gamma.star
	#(photorespiratory compensation point, Pa)
	if (coef.set == "Sharkey") {
	  T.Gamma.star = exp(11.187 - 24.46 / (R * TleafK))
	  }
	if (coef.set == "SilvaPerez") {
	  T.Gamma.star = (37.74 * exp(24.42 * (mean(Acidat[ , "Tleaf"]) - 25) /
	                                (R * 298 * TleafK))) * 1e-6 * 1e5
	  }
	                                  #microbar -> bar -> Pa
	
	if (any(is.null(Gamma.star), is.na(Gamma.star) & Gamma.star.fixed == TRUE)) {
	  fixed$Gamma.star = T.Gamma.star
	  }
	if (any(!is.na(Gamma.star) & Gamma.star.fixed == TRUE)) {
	  fixed$Gamma.star = Gamma.star
	  }
	if (any(is.na(Gamma.star) & Gamma.star.fixed != TRUE))	{
	  params$Gamma.star = T.Gamma.star
	  }
	if (any(!is.na(Gamma.star) & Gamma.star.fixed != TRUE)) {
	  params$Gamma.star = Gamma.star
	  }
	
	#Kco
	if (coef.set == "Sharkey") {
	  T.Ko = exp(12.3772 - 23.72 / (R * TleafK)) * 1e3 #oxygenation coefficient Pa
	  T.Kc = exp(35.9774 - 80.99 / (R * TleafK))#carboxylation coefficient Pa
	  }
								
	if (coef.set == "SilvaPerez") {
	  T.Ko = (166 * exp(33.6 * (mean(Acidat[ , "Tleaf"]) - 25) /
	                      (R * 298 * TleafK))) * 1e-3 * 1e5
	  #oxygenation coefficient mbar -> bar -> Pa
	  T.Kc = (272 * exp(93.72 * (mean(Acidat[ , "Tleaf"]) - 25) /
	                      (R * 298 * TleafK))) * 1e-6 * 1e5
	  #carboxylation coefficient microbar->bar->Pa
	  }
	
	T.Kco = T.Kc + T.Kc * (O2 / T.Ko)
	#print(T.Kco)
	
	if (any(is.null(Kco), is.na(Kco) & Kco.fixed == TRUE)) { fixed$Kco = T.Kco }
	if (any(!is.na(Kco) & Kco.fixed == TRUE)) { fixed$Kco = Kco }
	if (any(is.na(Kco) & Kco.fixed != TRUE))	{ params$Kco = T.Kco }
	if (any(!is.na(Kco) & Kco.fixed != TRUE)) { params$Kco = Kco }
	
	#Rd	
	## unused: to obtain Rd at Tleaf using Sharkey et al 2007
	## Rd.T<-Rd*exp(18.7145-46.39/(R*TleafK))
	#An Rd value *must* be supplied 
	if (any(is.null(Rd), is.na(Rd) & Rd.fixed != TRUE)) {params$Rd = 2}
	if (any(!is.na(Rd) & Rd.fixed != TRUE)) { params$Rd = Rd }
	if (any(is.na(Rd) & Rd.fixed == TRUE)) {
	  stop("Rd not provided and Rd.fixed set to TRUE")
	  }
	if (any(!is.na(Rd) & Rd.fixed == TRUE)) { fixed$Rd = Rd }
	
	#gm
	if (any(is.null(gm), is.na(gm) & gm.fixed != TRUE)) { params$gm = 5 }
	if (any(!is.na(gm) & gm.fixed == TRUE)) { fixed$gm = gm }
	if (all(is.na(gm)&gm.fixed==TRUE,var.J=F)) {
	  stop("neither of gm or var.J = TRUE provided and gm.fixed set to TRUE")
	  }
	if (any(!is.na(gm) & gm.fixed != TRUE)) { params$gm = gm }
	
	#via variable J - note that the equations here assume NADPH limitation
	if (var.J == TRUE){
	  if(!is.na(fixed$Rd) & !is.na(fixed$Gamma.star)){
	    gm.all <- Acidat$A /
	      (Acidat$Pci - (
	        (fixed$Gamma.star * (Acidat$ETR + 8 * (Acidat$A + fixed$Rd))) /
						                   (Acidat$ETR - 4 * (Acidat$A + fixed$Rd))
	        )
	       )
	    gm.accept <- 12 * fixed$Gamma.star * Acidat$ETR /
	      (Acidat$ETR - 4 * (Acidat$A + fixed$Rd))^2
	    gm.accept <- 1 < gm.accept & gm.accept < 5
	    #threshold from Harley et al, 1992. adjusted to Pa units
	    gm = mean(gm.all[gm.accept == TRUE], na.rm = TRUE)
	    if (any(is.na(gm), is.null(gm))) {
	      stop("variable J estimate of gm failed")
	      }
	    #Acidat$gm <- replace(gm.all * gm.accept,
			#                     (gm.all * gm.accept) == 0,
			#                     NA
			#                     )
	    fixed$gm = gm
	    } else {
	      stop("Rd & Gamma star not provided and var.J set to TRUE")
	    }
	  }
	
	#provide a dataframe and list of parameters for the fitting routine
	params = params
	if (!is.null(fixed)) {
	  fixed = as.list(unlist(fixed)[c(2:length(fixed))])
		list(data = Acidat,
		     params = params,
		     fixed = fixed,
		     stoich = stoich,
		     coef.set = coef.set,
		     physical = physical
		    )
	} else {
	    list(data = Acidat,
	         params = params,
	         fixed = NULL,
	         stoich = stoich,
	         coef.set = coef.set,
	         physical = physical
	    )
	  }
	}

#Function to generate and append a "limitation states" table
# input is output from Acifit.input
gen.lim.states <- function(input){
  #generates a matrix of potential limitation state combinations
  #with N depending on fixed versus estimated parameters.
	#Minimum N per segment
  #(discounting Rd and Gamma.star (not used in carboxylation state equations)):
	#if Kco is fixed:
	#	>=2 for Ac, >=2 Aj, and >=2 TPU
	#if Kco is to be estimated
	#	>=3 for Ac, >=2 Aj, and >=2 TPU
	#note: if the full TPU model is used, this will change
  # alpha will be an additional parameter for TPU
	
	#each row in the below is a 'state set'
	#generate all options (following Shaun Nielsen's lead)
	Nd <- nrow(input$data)
	a = rbind(rep(1, Nd), 1 * upper.tri(diag(x = Nd)))
	b = rbind(rep(2, Nd), 2 * upper.tri(diag(x = Nd)))
	mat12 = lapply(1:ncol(a),
	               function(.){ cbind(a[ , (1:.)], b[ , -(1:.)]) }
	               )
	grid.all <- do.call(rbind, mat12)
	grid.all <- unique(apply(grid.all, 1, paste,collapse = ''))
	grid.all <- do.call(rbind, strsplit(grid.all, ''))
	mode(grid.all) <- 'numeric'
	#The following hashed out expand.grid method failed when N >= 15 or so
	#grid.all <- expand.grid(split(rep(c(0:2),
	#                         each = nrow(input$data)),
	#                         seq(1, nrow(input$data), 1)),
	#                         KEEP.OUT.ATTRS = FALSE
	#                       )
	#use rules to narrow down to possible sets
	#1) ensure that states are always in correct, ascending, order
	#grid.all <- grid.all[!apply(
	#                             grid.all,
	#                             1,
	#                             function(.){
	#                               any(
	#                                    .[c(2:length(.))] -
	#                                    .[c(1:(length(.) - 1))] < 0)
	#                               }
	#                            ), ]
	#2) in this version, don't exclude sole limitation states
	#grid.all <- grid.all[!apply(
	#                             grid.all,
	#                             1,
	#                             function(.){
	#                               all(
	#                                   .[c(2:length(.))] -
	#                                   .[c(1:(length(.) - 1))] == 0)
	#                              }
	#                            ), ]
	#3) ensure correct N for fits
	if (is.null(input$params[["Kco"]])){
	  #3 a) Rubisco limited N = 0, or N > 1
	  grid.all <- grid.all[apply(grid.all,
	                             1,
	                             function(.){ length(.[. == 0]) != 1 }
	                             ), ]
	  } else {
	    #3 a) Rubisco limited N = 0, or N > 2
	    grid.all <- grid.all[apply(grid.all,
	                               1,
	                               function(.){ length(.[. == 0]) != 1 }
	                               ), ]
	    grid.all <- grid.all[apply(grid.all,
	                               1,
	                               function(.){ length(.[. == 0]) != 2}
	                               ), ]
	    }
	#3 b) RuBP limited 0 or n > 1
	grid.all <- grid.all[apply(grid.all,
	                           1,
	                           function(.){ length(.[. == 1]) != 1}
	                           ), ]
	#3 c) TPU limited 0 or n > 1
	grid.all <- grid.all[apply(grid.all,
	                           1,
	                           function(.){ length(.[. == 2]) != 1}
	                           ), ]
	
	row.names(grid.all) <- c(1:nrow(grid.all))
	
	input$lim.states <- grid.all
	input
	}

#Function that adapts total list of params from Acifit.input
#to those needed for a fit on a given line in the lim.states matrix
#also appends lim, a line from lim.states
#input is equivalent to output from Acifit.input
pick.params <- function(input, lim){
	
	lim <- lim
	pick.input <- input
	
	if (!any(lim == 0)) { 
						if (!is.null(pick.input$params[["Vcmax"]])) {
						  pick.input$params[["Vcmax"]] = NA
						  }
						if (!is.null(pick.input$params[["Kco"]])) {
						  pick.input$params[["Kco"]] = NA
						}
						if (!is.null(pick.input$fixed[["Vcmax"]])) {
						  pick.input$fixed[["Vcmax"]] = NA
						}
						if (!is.null(pick.input$fixed[["Kco"]])) {
						  pick.input$fixed[["Kco"]] = NA
						}
	}
	
	if (!any(lim == 1)) {
	  if (!is.null(pick.input$params[["J"]])) {
	    pick.input$params[["J"]] = NA
	  }
	  if (!is.null(pick.input$fixed[["J"]])) {
	    pick.input$fixed[["J"]] = NA
	  }
	}
	
	if (!any(lim == 2)) {
	  if (!is.null(pick.input$params[["TPU"]])) {
	    pick.input$params[["TPU"]] = NA
	  }
	  if (!is.null(pick.input$fixed[["TPU"]])) {
	    pick.input$fixed[["TPU"]] = NA
	  }
	}

	pick.input$params <- pick.input$params[!is.na(pick.input$params)]
	pick.input$fixed <- pick.input$fixed[!is.na(pick.input$fixed)]	
	
	list(data = pick.input$data,
	     params = pick.input$params,
	     fixed = pick.input$fixed,
	     stoich = pick.input$stoich,
	     lim = lim
	)
	
	}

#Function that evaluates Rubisco limited CO2 assimilation
#after Ethier & Livingstone/Gu
Rubisco.limited <- function(Vcmax, Kco, Gamma.star, Rd, gm, Pci.Rbc){
  b.Rbc = (Vcmax - Rd + (Pci.Rbc + Kco) * gm)
  c.Rbc = (((Pci.Rbc - Gamma.star) * Vcmax - (Pci.Rbc + Kco) * Rd) * gm)
  (b.Rbc - sqrt(b.Rbc^2 - 4 * c.Rbc)) / 2
}

#Function evaluating RuBP regeneration limited CO2 assimilation
#after Ethier & Livingstone/Gu
RuBP.limited <- function(J, Gamma.star, Rd, gm, Pci.J, stoich){
  if (stoich == "ATP") { c1 = 4.5; c2 = 10.5/4.5 } else { c1 = 4; c2 = 2 }
	b.J = ((J / c1) - Rd + ((Pci.J + (c2 * Gamma.star)) * gm))
	c.J = ((((Pci.J - Gamma.star) * J / c1) -
	          ((Pci.J + (c2 * Gamma.star)) * Rd)) * gm)
	(b.J - sqrt(b.J^2 - 4 * c.J)) / 2
}

#Function for TPU limited CO2 assimilation after Gu, with alpha = 0
TPU.limited <- function(TPU, Gamma.star, Rd, gm, Pci.TPU){
  b.P = (3 * TPU - Rd + (Pci.TPU - Gamma.star) * gm)
	c.P = (((Pci.TPU - Gamma.star) * 3 * TPU - (Pci.TPU - Gamma.star) * Rd) * gm)
	(b.P - sqrt(b.P^2 - 4 * c.P)) / 2
}

#Function generating cost function (residual ss)
#using parameters, data, a row from the lim state matrix, & stoichiometry
#appropriate input chosen from output of pick.params
#depends on functions: Rubisco.limited, RuBP.limited and TPU.limited
#modified 0819 to enable co-limited points to be fit
Aci.cost <- function(params, fixed, data, lim, stoich){
  #lim is the single row within the lim.states to which model is fit
	#par has to be entered separately so optim works
	#inputs must be as named vectors so that optim works
								
	#joint parameters
	Rd = ifelse(!is.na(params["Rd"]),
	            params["Rd"],
	            fixed["Rd"]
	            )
	Gamma.star = ifelse(!is.na(params["Gamma.star"]),
	                    params["Gamma.star"],
	                    fixed["Gamma.star"]
	                    )
	gm = ifelse(!is.na(params["gm"]),
	            params["gm"],
	            fixed["gm"]
	            )
	
	#Rubisco limited
	Rubisco.resid = NA
	if (!is.na(params["Vcmax"])) {
	  Vcmax = params["Vcmax"]
	  Kco = ifelse(!is.na(params["Kco"]),
	               params["Kco"],
	               fixed["Kco"]
	               )
	  Pci.Rbc = data[lim == 0 | lim == 3 | lim == 5, "Pci"]
	  
	  Rubisco.fitted = Rubisco.limited(Vcmax, Kco, Gamma.star, Rd, gm, Pci.Rbc)
	  Rubisco.resid = data[lim == 0 | lim == 3 | lim == 5, "A"] - Rubisco.fitted
	  }
	
	#RuBP regeneration limited
	RuBP.resid = NA
	if (!is.na(params["J"])) {
	  J = params["J"]
	  Pci.J <- data[lim == 1 | lim == 3 | lim == 4, "Pci"]
	  
	  RuBP.fitted = RuBP.limited(J, Gamma.star, Rd, gm, Pci.J, stoich)
	  RuBP.resid = data[lim == 1 | lim == 3 | lim == 4, "A"] - RuBP.fitted
	  }
	
	#TPU limited
	TP.resid = NA
	if (!is.na(params["TPU"])) {
	  TPU = params["TPU"]
	  Pci.TPU = data[lim == 2 | lim == 4 | lim == 5, "Pci"]
	  
	  TP.fitted = TPU.limited(TPU, Gamma.star, Rd, gm, Pci.TPU)
	  TP.resid = data[lim == 2 | lim == 4 | lim == 5, "A"] - TP.fitted
	  }
	
	#cost function					
	1 / 2 * sum(Rubisco.resid^2, RuBP.resid^2, TP.resid^2, na.rm = TRUE)
	#note: points allocated to two limitation states are double counted
	#so co-limited models are penalised!
	
}

#Function for nonlinear fit optimisation
#0819 - re-wrote to allow this to be slotted into
# either initial or swinging point (co-limited) routines
#input is the usual object, p.para is pick.params output
opt.pars <- function(input, p.para){
  
  p.p = unlist(p.para$params)
	f.p = unlist(p.para$fixed)
	d.p = p.para$data
	l.p = p.para$lim
	upp = rep(Inf, length(p.para$params))
	low = rep(0, length(p.para$params))
	stoich = input$stoich
	
	#obtain an unbounded estimate using a non-derivative based method
	#optim(p.p,
	#		fn = Aci.cost,
	#		method = "Nelder-Mead",
	#		#upper = upp, lower = low,
	#		control = list(maxit = 100000),
	#		fixed = f.p, data = d.p, lim = l.p, stoich = stoich
	#		)
	
	optim(p.p,
	      fn = Aci.cost,
	      fixed = f.p,
	      data = d.p,
	      lim = l.p,
	      stoich = stoich,
	      method = "L-BFGS-B",
	      upper = upp, lower = low,
	      control = list(maxit = 10000)
	)
	
	}

#Function taking Acifit.input pre-processed by gen.lim.states
#and using opt.pars to implement pick.params and Aci.cost,
#providing fitted models for every line in lim.states.
#Appends estimated parameters to input
Aci.cost.fits <- function(input){
  
  lim.states <- input$lim.states
  stoich = input$stoich
	
	#dataframe to hold output for each individual fit
	cost <- data.frame(matrix(NA, nrow = nrow(lim.states), ncol = 9))
	names(cost) <- c("convergence", "cost", "Vcmax", "J", "TPU", "Rd",
	                 "Gamma.star", "Kco", "gm"
	                 )
	
	#loop fitting models and calculating output for all rows in lim.states
	for (k in 1:nrow(lim.states)){
	  print(k)
	  p.para = pick.params(input, lim.states[k, ])
	  fit <- opt.pars(input, p.para)
	  
	  cost[k, "convergence"] = fit$convergence
	  cost[k, "cost"] = fit$value
	  
	  for (l in names(p.para$params)){
	    cost[k, l] <- fit$par[l]
	    }
	  
	  for (m in names(p.para$fixed)){
	    cost[k, m] <- p.para$fixed[m]
	    }
	  
	  }
	
	input$models = cost
	input
	
	}

#function updating input$models (i.e., output from Aci.cost.fits)
#by adding a column indicating whether models are admissible
#based on whether minimum carboxylation rate (depending on Cc)
#fits with assigned limitation state (depending on Ci)
#appends the result to test.cost,
#and adds the lim.states and pred.states to data
Acifit.admit <- function(input){
  
  cost <- input$models
  lim.states <- input$lim.states
  pred.states <- data.frame(matrix(NA, nrow(lim.states), ncol(lim.states)))
  stoich = input$stoich
  cost$admit <- rep(NA, nrow(cost))
  
  attach(cost)
  
  for (i in 1:nrow(cost)){
    #net assimilation preditions for the entire range of Pci in the data
		Ac = Rubisco.limited(Vcmax[i],
		                     Kco[i],
		                     Gamma.star[i],
		                     Rd[i],
		                     gm[i],
		                     input[["data"]]$Pci
		                     )
		
		Aj = RuBP.limited(J[i],
		                  Gamma.star[i],
		                  Rd[i],
		                  gm[i],
		                  input[["data"]]$Pci,
		                  stoich
		                  )
		
		Ap = TPU.limited(TPU[i],
		                 Gamma.star[i],
		                 Rd[i],
		                 gm[i],
		                 input[["data"]]$Pci
		                 )
		
		#calculate Cc based on A and gm
		Cc = input[["data"]]$Pci - Ac / gm[i]
		Cj = input[["data"]]$Pci - Aj / gm[i]
		Cp = input[["data"]]$Pci - Ap / gm[i]
		
		#carboxylation rates
		Wc = Vcmax[i] * Cc / (Cc + Kco[i])
		if (stoich == "ATP") { c1 = 4.5; c2 = 10.5 } else { c1 = 4; c2 = 8 }
		Wj = J[i] * Cj / (c1 * Cj + c2 * Gamma.star[i])
		Wp = 3 * TPU[i] * Cp / (Cp - Gamma.star[i])
		
		#set values for TPU that fall below threshold Cc:
		#(Cc>(1+3*alpha)*Gamma.star) so that they will be safely excluded.
		#Rounding to a sensible accuracy is necessary here to select correctly
		Wp <- replace(Wp, round(Cp, 4) <= round(Gamma.star[i], 4), Inf)
		
		#same, just in case there are any NaNs generated for any of the states
		#0121 noticed that these lines were nonsense...
		#nonsense lines are shown here at right, hashed out
		Wp <- replace(Wp, is.na(Wp), Inf) #Wp <- replace(Wp, Wp == NA, Inf)
		Wc <- replace(Wc, is.na(Wc), Inf) #Wc <- replace(Wc, Wc == NA, Inf)
		Wj <- replace(Wj, is.na(Wj), Inf) #Wj <- replace(Wj, Wj == NA, Inf)
		
		#all present
		if (all(!is.na(Wc)) & all(!is.na(Wj)) & all(!is.na(Wp))) {
		  pred.states[i, ] <- replace(pred.states[i, ], Wc < Wj & Wc < Wp, 0)
		  pred.states[i, ] <- replace(pred.states[i, ], Wj < Wp & Wj < Wc, 1)
		  pred.states[i, ] <- replace(pred.states[i, ], Wp < Wc & Wp < Wj, 2)
		  }
		
		#only Wc
		if (all(!is.na(Wc)) & all(is.na(Wj)) & all(is.na(Wp) | Wp == Inf)) {
		  pred.states[i, ] <- rep(0, ncol(pred.states))
		  }
		
		#only Wj
		if (all(!is.na(Wj)) & all(is.na(Wc)) & all(is.na(Wp) | Wp == Inf)) {
		  pred.states[i, ] <- rep(1, ncol(pred.states))
		  }
		
		#only Wp
		if (all(!is.na(Wp)) & all(is.na(Wc)) & all(is.na(Wj))) {
		  pred.states[i, ] <- replace(pred.states[i, ], Wp != Inf, 2)
		  pred.states[i, ] <- replace(pred.states[i, ], Wp == Inf, Inf)
		  }
		
		#Wc & Wj
		if (all(!is.na(Wc)) & all(!is.na(Wj)) & all(is.na(Wp) | Wp == Inf)) {
		  pred.states[i, ] <- replace(pred.states[i, ], Wc < Wj, 0)
		  pred.states[i, ] <- replace(pred.states[i, ], Wj < Wc, 1)
		  }
		
		#Wc & Wp
		if (all(!is.na(Wc)) & all(!is.na(Wp)) & all(is.na(Wj))) {
		  pred.states[i, ] <- replace(pred.states[i, ], Wc < Wp, 0)
		  pred.states[i, ] <- replace(pred.states[i, ], Wp < Wc, 2)
		  }
		
		#Wj & Wp
		if (all(!is.na(Wj)) & all(!is.na(Wp)) & all(is.na(Wc))) {
		  pred.states[i, ] <- replace(pred.states[i, ], Wj < Wp, 1)
		  pred.states[i, ] <- replace(pred.states[i, ], Wp < Wj, 2)
		  }
		
		lim.sti <- as.numeric(lim.states[i, ])
		pred.sti <- as.numeric(pred.states[i, ])
		
		cost[i, "admit"] <- all(lim.sti == pred.sti)
		
		#0819 added the following to ensure that failures of pred.states
		# are counted as non-admissible
		#This will happen if, e.g. -ve Vcmax, J, or TPU
		#it therefore isn't an issue if L-BFGS-B is used
		#and should not affect this code unless op.pars is modified
		cost[i, "admit"] <- replace(cost[i, "admit"],
		                            is.na(cost[i, "admit"]),
		                            FALSE
		                            )
		
		}
	
	input$models = cost
	input$pred.states = pred.states
	detach(cost)
	input

	}
	
#Function added 0819 to handle 'swinging-points' identification sensu Gu et al
#- finds rows where lim.states and pred.states differ = inadmissibile
#- finds models with initial limitation states matching predicted states
#   from inadmissible models
#- checks that initial and predicted are reciprocal in this pairing (swinging)
#- if limitation state is swinging for one point:
#   - a row is added to the limitation state matrix
#      with the average limitation state assigned to the swinging point
#				(this allows easy boolean assignment of the states later)
#				the model is refit with the point allocated to both states at once
#				admissibility of the swinging point model is set to TRUE
#- if the above process does not work, model admissibility remains FALSE
#this last point is different from the wording of Gu et al.,
#but they seem to indicate that this is effectively what they are doing
Acifit.colimit <- function(input){
  
  models <- input$models
	lim.states <- data.frame(input$lim.states)
	pred.states <- data.frame(input$pred.states)
	stoich = input$stoich
	inadm =! input$models$admit
	na.pred = apply(!is.na(pred.states), 1, any)
	
	to.check <- lim.states[inadm & na.pred, ]
	
	for (i in 1:nrow(to.check)){
	  
	  pred.inad = pred.states[as.numeric(rownames(to.check)[i]), ]
		lim.alt.id = apply(lim.states, 1, function(.){ all(. == pred.inad) })
		lim.alt = lim.states[lim.alt.id, ]
		pred.alt = pred.states[lim.alt.id, ]
		
		if (!nrow(pred.inad) == 0 &
		    !all(is.na(pred.inad)) &
		    !nrow(lim.alt) == 0 &
		    !all(is.na(lim.alt)) &
		    !nrow(pred.alt) == 0 &
		    !all(is.na(pred.alt))
		    ) {
		  
		  #081913 updated logic gating below to prevent duplication
		  #(since we know anything passing the first test has a mirror)
			if (all(to.check[i, ] == pred.alt)) {
			  
			  #indexing for new rows to be added to lim.states, pred.states, & models
				new.row = nrow(lim.states) + 1
				
				#generate a colimited state vector
				colims <- rep(NA, length(to.check[i, ]))
				for (j in 1:length(colims)){
					if (to.check[i, j] == lim.alt[j]) {
					  colims[j] <- to.check[i, j]
					  } else {
					    
					    if ((to.check[i, j] == 0 & lim.alt[j] == 1) |
					        (to.check[i, j] == 1 & lim.alt[j] == 0)
					    ) { colims[j] <- 3 }
					    
					    if ((to.check[i, j] == 1 & lim.alt[j] == 2) |
					      (to.check[i, j] == 2 & lim.alt[j] == 1)
					    ) { colims[j] <- 4 }
				    
					    if ((to.check[i, j] == 0 & lim.alt[j] == 2) |
				        (to.check[i, j] == 2 & lim.alt[j] == 0)
					    ) { colims[j] <- 5 }
					    
					  }
				  }
				
				#before extending the tables,
				# check that this isn't duplicating a row that has already been added
				if (!all(lim.states[nrow(lim.states), ] == colims)){
					lim.states[new.row, ] <- colims
					pred.states[new.row, ] <- colims
					models[new.row, ] <- rep(NA, ncol(models))
			
					#refit the model
					p.para <- pick.params(input, colims)
					
					fit <- opt.pars(input, p.para)
					
					models[new.row, "convergence"] = fit$convergence
					models[new.row, "cost"] = fit$value
										
					for (k in names(p.para$params)){
						models[new.row, k] = fit$par[k]
						}
					
					for (l in names(p.para$fixed)){
						models[new.row, l] = p.para$fixed[l]
						}

					#and assign it as admissible
					# (assured by the definition of the swinging point)	
					models[new.row, "admit"] = TRUE
				}
			}
		}
		}
	
	#return values as part of the input object
	input$models = models
	input$lim.states = lim.states
	input$pred.states = pred.states
	input
	
}


#functions for additional parameters to be derived from individual models
 
#Compensation point

#assuming Rubisco limited carboxylation 
#sign in numerator corrected compared with earlier scripts
A.0.Ac <- function(Vcmax, Gamma.star, Rd, Kco){
  (Vcmax * Gamma.star + Rd * Kco) / (Vcmax - Rd)
  }

#assuming RuBP regeneration limited carboxylation 
A.0.Aj <- function(J, Gamma.star, Rd, stoich){
  if (stoich == "ATP") { c1 = 4.5; c2 = 10.5 } else { c1 = 4; c2 = 8 }
  (Gamma.star * (J + c2 * Rd)) / (J - c1 * Rd)
}

#TPU cannot be limiting at the compensation point because
# Wp only applies if Cc > (1 + 3 * alpha) * Gamma.star
# when alpha = 0, the compensation point for TPU == Gamma.star
# if alpha > 0, it is < Gamma.star

#Transition Ci for Ac == Aj
# from quadratic equality derived from Cc level Ac and Aj
Ac.Aj <- function(Vcmax, J, Gamma.star, Kco, gm, Rd, stoich){
  if (stoich == "ATP") { c1 = 4.5; c2 = 10.5 } else { c1 = 4; c2 = 8 }
  a = c1 * Vcmax - J
  b = (J * Gamma.star) + ((c2 - c1) * Vcmax * Gamma.star) - (J * Kco)
  c = (J * Kco * Gamma.star) - (c2 * Vcmax * Gamma.star^2)
  roots <- c((-b - sqrt(b^2 - 4 * a * c)) / (2 * a),
             (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
             )
  Pcc = max(roots)
  Ac.Aj = Vcmax * ((Pcc - Gamma.star) / (Pcc + Kco)) - Rd #Ac predicted from Cc
	Pcc + (Ac.Aj / gm)
}

#transition ci for Aj=Ap from quadratic equality derived from Cc level Aj and Ap 
Aj.Ap <- function(J, TPU, Gamma.star, gm, Rd, stoich){ #alpha = 0
  if (stoich == "ATP") { c1 = 4.5; c2 = 10.5 } else { c1 = 4; c2 = 8 }
  a = J - c1 * 3 * TPU
  b = Gamma.star * (3 * TPU * (c1 - c2) - 2 * J)
   #if alpha != 0: Gamma.star * (3 * TPU * (c1 - c2) - J * (2 + 3 * alpha))
  c = Gamma.star^2 * (J + c2 * 3 * TPU)
   #if alpha !=0: Gamma.star^2 * (J * (1 + 3 * alpha) + c2 * 3 * TPU)
	roots <- c((-b - sqrt(b^2 - 4 * a * c)) / (2 * a),
	           (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
	           )
	Pcc = max(roots)
	Aj.Ap = J * ((Pcc - Gamma.star) / (c1 * Pcc + c2 * Gamma.star)) - Rd
	 #Aj predicted from Cc
	Pcc + (Aj.Ap / gm)
}

#Transition Ci for Ac == Ap from quadratic equality of Cc level Ac and Ap 
Ac.Ap <- function(Vcmax, Kco, TPU, Gamma.star, gm, Rd){ #alpha = 0
  a = 3 * TPU - Vcmax
  b = 2 * Vcmax * Gamma.star + 3 * TPU * (Kco - Gamma.star)
   #if alpha != 0:
   #Vcmax * Gamma.star * (2 + 3 * alpha) + 3 * TPU * (Kco - Gamma.star)
	c = -Gamma.star * (Vcmax * Gamma.star + 3 * TPU * Kco)
	 #if alpha != 0:
	 #-Gamma.star * (Vcmax * Gamma.star * (1 + 3 * alpha) + 3 * TPU * Kco)
	roots <- c((-b - sqrt(b^2 - 4 * a * c)) / (2 * a),
	           (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
	           )
	Pcc = max(roots)
	Ac.Ap = Vcmax * ((Pcc - Gamma.star) / (Pcc + Kco)) - Rd #Ac predicted from Cc
	print(Ac.Ap)
	Pcc + (Ac.Ap / gm)
}

#There's no way based on Sharkey/Bernacchi to normalise everything back to 25 C,
#just the primary coefficients

#Vcmax at 25 C after Sharkey/Bernacchi
Vcmax.25 <- function(Vcmax, R, Tleaf.K){
  Vcmax * exp(26.355 - 65.33 / (R * Tleaf.K))
  }

#J at 25 C after Sharkey/Bernacchi
#(need to investigate which stoichiometry is relevant here.
# The Sharkey 20017 paper will probably use 4/8,
# but the original derivation in Bernacchi might be 4.5/10.5)
J.25 <- function(J, R, Tleaf.K){
  J * exp(17.71 - 43.9 / (R * Tleaf.K))
  }

#TPU at 25 C after Sharkey/Bernacchi
#(originally Harley 1992, and applying form of equation used for gm below)
TPU.25 <- function(TPU, R, Tleaf.K){
  TPU * (
    exp(21.46 - 53.1 / (R * Tleaf.K)) /
           (1 + exp((0.65 * Tleaf.K - 201.8) / (R * Tleaf.K)))
    )
}

#Rd at 25 C after Sharkey/Bernacchi
Rd.25 <- function(Rd, R, Tleaf.K, coef.set = "Sharkey"){
  
  if (coef.set == "Sharkey") {
    Rd.25 <- Rd * exp(18.7145 - 46.39 / (R * Tleaf.K))
  }
  
  if (coef.set == "SilvaPerez") {
    Rd.25 <- Rd / exp(46.39 * (Tleaf.K - 273.15 - 25) / (R * 298 * Tleaf.K))
    }
  
  Rd.25
  }

#gm at 25 C after Bernacchi (Sharkey equation is wrong)
gm.25 <- function(gm, R, Tleaf.K){
  gm * (
    exp(20.01 - 49.6 / (R * Tleaf.K)) /
      (1 + exp((1.4 * Tleaf.K - 437.4) / (R * Tleaf.K)))
  )
}

#Using input list object with with models,
#and a fitted model selected from input$models
# (output of Aci.cost.fits with or without admit added)
#- Calculate parameters
#   including temperature corrections.
#- Prepare a complete set of information for the model
#   that can be used in a plot function
#- Optionally, if steady state values for A and Pci are available
#   these can be added to calculate Farquhar&Sharkey stomatal limitation
#   - corrected this 0320 to work in situations with no A.op or Pci.op
#   - further corrected 0720 to allow SS ca (Pca.op) to be specified
#       to give a more correct Farquhar&Sharkey calc
Acifit.evaluate <- function(model,
                            input,
                            A.op = NULL,
                            Pci.op = NULL,
                            Pca.op = NULL
                            ){
  stoich = input$stoich
	attach(model)
	
	#Ac, Aj, Ap over physiologically relevant Cc
	Pci <- round(seq(0, 150, 0.1), 1) #round necessary to prevent failures in A.ca
	Ac <- Rubisco.limited(Vcmax, Kco, Gamma.star, Rd, gm, Pci)
	Aj <- RuBP.limited(J, Gamma.star, Rd, gm, Pci, stoich)
	Ap <- TPU.limited(TPU, Gamma.star, Rd, gm, Pci)
	
	#Assimilation rate at Ci=40 Pa
	A.ca <- min(Ac[Pci == round(Pca.op, 1)],
	            Aj[Pci == round(Pca.op, 1)],
	            Ap[Pci == round(Pca.op, 1)],
	            na.rm = TRUE
	            )
	
	#Farquhar and Sharkey stomatal limitation
	if (!is.null(A.op)){
	  Ls = (A.ca - A.op) / A.ca
	  }
	
	#compensation point
	if (!is.na(Vcmax)) {
	  A.0 = A.0.Ac(Vcmax, Gamma.star, Rd, Kco)
	  } else {
	    if (!is.na(J)) {
	      A.0 = A.0.Aj(J, Gamma.star, Rd, stoich)
	    } else {
	      A.0 = NA
	    }
	    }
	
	#Ac/Aj transition ci
	if (all(!is.na(Vcmax), !is.na(J))) {
	  Ac.Aj = Ac.Aj(Vcmax, J, Gamma.star, Kco, gm, Rd, stoich)
	} else {
	    Ac.Aj = NA
	    }
	
	#Aj/Ap transition ci
	if (all(!is.na(J), !is.na(TPU))) {
	  Aj.Ap = Aj.Ap(J, TPU, Gamma.star, gm, Rd, stoich)
	} else {
	    Aj.Ap = NA
	    }
	
	#Ac/Ap transition ci
	if (all(!is.na(Vcmax), !is.na(TPU), is.na(J))) {
	  Ac.Ap = Ac.Ap(Vcmax, Kco, TPU, Gamma.star, gm, Rd)
	} else {
	    Ac.Ap = NA
	    }
	
	#temperature corrected values
	TleafK = input[["physical"]]$Tleaf.K
	R = input[["physical"]]$R
	coef.set = input[["coef.set"]]
	
	Vcmax.25 = Vcmax.25(Vcmax, R, TleafK)
	J.25 = J.25(J, R, TleafK)
	TPU.25 = TPU.25(TPU, R, TleafK)
	Rd.25 = Rd.25(Rd, R, TleafK, coef.set)
	gm.25 = gm.25(gm, R, TleafK)
	
	match.lims <- t(apply(input$models, 1, function(.){ . == model }))
	match.which <- apply(match.lims, 1, function(.){ all(is.na(.) | T == .) })
	chosen.lims = unlist(input$lim.states[match.which, ])
	names(chosen.lims) <- NULL
	
	input$chosen.lims = chosen.lims
	input$chosen.mod = model
	input$predicts = data.frame(Pci = Pci, Ac = Ac, Aj = Aj, Ap = Ap)
	input$F_Slimitation = c(A.ca = A.ca,
	                        A.op = A.op,
	                        Ls = Ls,
	                        Pci.op = Pci.op,
	                        Pca.op = Pca.op
	                        )
	input$ci.transitions = c(Gamma = A.0,
	                         Ac.Aj = Ac.Aj,
	                         Aj.Ap = Aj.Ap,
	                         Ac.Ap = Ac.Ap
	                         )
	input$temperature.normalised = c(Vcmax.25 = Vcmax.25,
	                                 J.25 = J.25,
	                                 TPU.25 = TPU.25,
	                                 Rd.25 = Rd.25,
	                                 gm.25 = gm.25
	                                 )
	
	detach(model)
	input
}


#plot the data and fitted model,
#identifying limitation states, and, if present, co-plotting phi[PSII]
#uses ouput from Acifit.evaluate
# i.e., all input, with parameters derived from a chosen model
#Adapted to insert the swinging point as co-limited green
plot.Acifit <- function(input,
                        xlim = NULL,
                        ylim = NULL,
                        PhiPSII = FALSE,
                        F_Slimitation = FALSE,
                        add.legend = TRUE,
                        ...
                        ){
  
  par(mar = c(5.5, 5.5, 3, 1), las = 1)
  
  lims.pal <- data.frame(lims = c(0:5),
                         cols = c(gray(1), gray(0), gray(0.7), rep("green", 3))
                         )
  
  col.picker <- function(.){
    for (i in 1:length(.)){
      .[i] <- as.character(lims.pal[lims.pal$lims == .[i], "cols"])
      }
    .
    }
	
	if (is.null(ylim)) {
	  y.dat <- input$data[, "A"]
	  ylim = range(y.dat[is.finite(y.dat)])
	  }
  
  if (is.null(xlim)) {
    x.dat <- input$data[, "Pci"]
    xlim = range(x.dat[is.finite(x.dat)])
    }
  
  if (PhiPSII == TRUE) {
    par(mar = c(5.5, 5.5, 3, 5.5), las = 1)
    }
  
  plot(A ~ Pci,
       data = input$data,
       xlim = xlim,
       ylim = ylim,
       type = "n",
       xlab = expression(italic(c)[i]~~(Pa)),
       ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1))
       )
  abline(0, 0)
  lines(Ac ~ Pci, data = input$predicts, lty = 2)
	lines(Aj ~ Pci, data = input$predicts, lty = 5)
	lines(Ap ~ Pci, data = input$predicts, lty = 6)
	lines(rep(input$ci.transitions["Gamma"], 2), c(-1e6, 1e6), lty = 3)
	lines(rep(input$ci.transitions["Ac.Aj"], 2), c(-1e6, 1e6), lty = 3)
	lines(rep(input$ci.transitions["Aj.Ap"], 2), c(-1e6, 1e6), lty = 3)
	lines(rep(input$ci.transitions["Ac.Ap"], 2), c(-1e6, 1e6), lty = 3)
	
	points(A ~ Pci,
	       data = input$data,
	       pch = 21,
	       bg = col.picker(input$chosen.lims)
	       )
	
	if (F_Slimitation == TRUE) {
	  points(input$F_Slimitation["Pci.op"], input$F_Slimitation["A.op"],
	         pch = 21,
	         cex = 1.2,
	         col = 4,
	         bg = 4
	         )
	  
	  points(input$F_Slimitation["Pca.op"], input$F_Slimitation["A.ca"],
	         pch = 21,
	         cex = 1.2,
	         col = 2,
	         bg=4
	         )
	  
	  lines(rep(input$F_Slimitation["Pca.op"], 2),
	        c(0, input$F_Slimitation["A.ca"]),
	        col = 4
	        )
	  
	  lines(rep(input$F_Slimitation["Pca.op"] - 0.5, 2),
	        c(input$F_Slimitation["A.op"], input$F_Slimitation["A.ca"]),
	        col = 2,
	        lwd = 2
	  )
	}
	
	if (PhiPSII == TRUE) {
	  phi.ex <- max(ylim) / (1.2 * max(input$data[ , "PhiPS2"]))
	  
	  points(phi.ex * PhiPS2 ~ Pci,
	         data = input$data,
	         pch = 24,
	         bg = col.picker(input$chosen.lims)
	         )
	  
	  axis(side = 4,
	       at = seq(round(min(ylim), -1), round(max(ylim), -1), 10),
	       labels = round(
	         seq(round(min(ylim), -1),
	             round(max(ylim), -1), 10) / phi.ex,
	         1
	         )
	       )
	  
	  mtext(side = 4, expression(phi[PSII]), las = 3, line = 3, cex = 1)
	  
	  if (add.legend == TRUE) {
	    legend(0.95 * diff(xlim) + min(xlim), 0.05 * diff(ylim) + min(ylim),
	           xjust = 1,
	           yjust = 0,
	           bg = gray(1),
	           pch = c(21, 24, rep(22, 3)),
	           pt.bg = c(rep(0, 2), col.picker(c(0:2))),
	           lty = c(rep(NA, 2), 2, 5, 6),
	           legend = expression(italic(A),
	                               italic(Phi)[PSII],
	                               italic(A)[C],
	                               italic(A)[J],
	                               italic(A)[P]
	                               )
	    )
	    }
	  
	  } else {
	    if (add.legend == TRUE) {
	      legend(0.95 * diff(xlim) + min(xlim),
	             0.05 * diff(ylim) + min(ylim),
	             xjust = 1,
	             yjust = 0,
	             bg = gray(1),
	             pch = rep(21, 3),
	             pt.bg = col.picker(c(0:2)),
	             lty = c(2, 5, 6),
	             legend = expression(italic(A)[C],
	                                 italic(A)[J],
	                                 italic(A)[P]
	                                 )
	      )
	    }
	  }
}
