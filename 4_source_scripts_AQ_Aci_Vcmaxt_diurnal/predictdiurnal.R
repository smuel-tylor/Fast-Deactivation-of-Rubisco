#Model of photosynthetic response to PPFD
#follows 031907predictdiurnal.R

#0420 load AQfit for convenience function
#0121 corrected this to test for function first
# and ensure most up to date version of AQfit is being used!
if (!exists("AQ.form", mode = "function")) {
  source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/AQfit.R"))
}
  
#Function to estimate tau.down
#Af = final A from induction fit
#Ai = initial A from induction fit
#shade.s = duration of imposed shade in seconds
#0221 modified as part of revisions to reflect decline to Ai, rather than 0
tau.down <- function(Af, Ai, shade.s){
  shade.s / log(-(Ai - Af))
}

#Diurnal model used with Zhu et al 2006 PPFD simulations
#light = diurnal light regime
#plant = light response and tau parameters as named vector
diurnal.mod <- function(light, plant){
  
  level2 <- light
  
  #parameters
	phi = plant["phi"]
	Asat = plant["Asat"]
	theta = plant["theta"]
	Rd = plant["Rd"]
	tau.up = plant["tau.up"]
	tau.down = plant["tau.down"]

	#expected net CO2 assimilation if this responds instantaneously to PPFD
	#(used gross assimilation A + Rd in Taylor & Long 2017
	# because Woodrow/Mott equation used to predict kinetics of induction
	# (A* = A * (Ci - Gamma)/(Ci* - Gamma) OR A* = (A + Rd) * Ci / Ci*
	#  i.e. assuming a linear response of A through the A/Ci intercept)
	# here, my tau etc. are from one point A/ci
	# so net CO2 assimilation is appropriate)
	level2$Af <- AQ.form(phi, Asat, theta, Rd, level2$PPFD)
	head(level2$Af)

	#timestep change in PPFD
	level2$dPPFD <- c(NA,
	                  level2$PPFD[2:nrow(level2)] -
	                    level2$PPFD[1:(nrow(level2) - 1)]
	                  )
	
	#identify if stepping up or not
	level2$PPFD.TF <- level2$PPFD > 0
	level2$dPPFD.TF <- level2$dPPFD != 0

	#duration of timestep and rules for carrying out operations
	#note, the timings here are specific to the Zhu simulation
	level2$t.since.dPPFD <- rep(NA, nrow(level2))
	for (i in 2:nrow(level2)){
	  nest.t <- ifelse(level2$PPFD.TF[i] & level2$dPPFD.TF[i],
	                   57.6,
	                   level2$dPPFD.TF[i - 1] + 57.6
	                   )
	  level2$t.since.dPPFD[i] <- ifelse(level2$PPFD.TF[i],
	                                    nest.t,
																			0
	  )
	}
	
	level2$t.since.dPPFD.TF <- (level2$t.since.dPPFD != 0 &
	                              !is.na(level2$t.since.dPPFD)
	)
	
	tf.1 <- level2$t.since.dPPFD.TF[c(1:(nrow(level2) - 1))]
	tf.2 <- level2$t.since.dPPFD.TF[c(2:nrow(level2))]
	level2$t.since.dPPFD.notmax <- c(NA, tf.1 < tf.2)
	
	level2$t.since.dPPFD.mid <- (level2$t.since.dPPFD.notmax &
	                               level2$t.since.dPPFD.TF)
	
	level2$dur.interval.s <- ifelse(level2$t.since.dPPFD.mid,
	                                0,
	                                level2$t.since.dPPFD
	                                )
	
	level2$do.calc <- c("TRUE", level2[c(2:nrow(level2)), "dur.interval.s"] > 0)

	#initial state for start of each timestep
	# Both up and down possibilities are predicted from previous timestep,
	# then correct direction (up or down) selected
	# based on Af versus Ai for the preceding step
	level2$Ai.up <- c(rep(0, nrow(level2)))
	level2$Ai.down <- c(rep(0, nrow(level2)))
	level2$Ai <- c(rep(0, nrow(level2)))
	
	for (i in 2:nrow(level2)){

		L2.Afp = level2[i - 1, "Af"]
		L2.t = i - (level2[i, "t.since.dPPFD"] / 57.6)

		level2[i, "Ai.up"] <- (L2.Afp -
		                         (L2.Afp - level2[L2.t, "Ai"]) *
		                         exp(-level2[L2.t, "t.since.dPPFD"] / tau.up)
		)
		
		level2[i, "Ai.down"] <- (L2.Afp -
		                           (L2.Afp - level2[L2.t, "Ai"]) *
		                           exp(-level2[L2.t, "t.since.dPPFD"] / tau.down)
		)
		
		level2[i, "Ai"] <- ifelse(level2[L2.t, "Af"] > level2[L2.t, "Ai"],
															level2[i, "Ai.up"],
															level2[i, "Ai.down"]
		)
		
	}
	
	#establish which timesteps show a difference between
	# instantaneous and rate-limited kinetics
	level2$Af.Ai <- rep(0, nrow(level2))

	for (i in 1:nrow(level2)){
		if (level2[i, "do.calc"] == TRUE) {
		  level2[i, "Af.Ai"] <- level2[i, "Af"] - level2[i, "Ai"]
		}
	  }
	
	#integrated gross assimilation

	#using rubisco limitation: applies only when PPFD is increasing
	level2$Ahat <- (level2$Af *
	                  level2$t.since.dPPFD - level2$Af.Ai *
	                  tau.up + level2$Af.Ai *
	                  tau.up *
	                  exp(-level2$t.since.dPPFD / tau.up)
	                )
	
	#Assimilation tracking PPFD
	level2$Aft <- level2$Af * level2$t.since.dPPFD

	level2$predict <- rep(0, nrow(level2))

	for (i in 2:nrow(level2)){
		level2[i, "predict"] <- ifelse(level2[i, "Ahat"] < level2[i, "Aft"],
		                               level2[i, "Ahat"],
		                               level2[i, "Aft"]
		)
		}
	
	#back converting to get instantaneous rates
	level2$F <- level2$Aft - level2$predict

	level2$predict.A <- level2$predict / level2$t.since.dPPFD
	level2$F.A <- level2$F / level2$t.since.dPPFD

	#output
	list(model = level2, params = plant)
}
