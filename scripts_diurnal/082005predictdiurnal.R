#Model of photosynthetic response to PPFD
#follows 031907predictdiurnal.R

#load AQfit for convenience function - only modification made 30 Apr '20
#source("C:/Users/taylor53/Box Sync/RIPE_Lancaster/Cowpea_GasExchange/Rscripts/011921AQfit.R")
source("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rscripts/011921AQfit.R")

#function to estimate tau.down
#Af = final A from induction fit
#Ai = initial A from induction fit
#shade.s = duration of imposed shade in seconds
tau.down <- function(Af, Ai, shade.s){
	-shade.s / log(Ai / Af)
}

diurnal.mod<-function(light, plant){#light is the diurnal light regime and plant the light response and tau parameters as named vector (see below)
	
	level2 <- light

	#parameters
	phi = plant["phi"]
	Asat = plant["Asat"]
	theta = plant["theta"]
	Rd = plant["Rd"]
	tau.up = plant["tau.up"]
	tau.down = plant["tau.down"]

	#expected A with instantaneous kinetics (was A+Rd for Woodrow/Mott, but my tau etc. are from A/ci, therefore Rd corrected)
	level2$Af <- AQ.form(phi, Asat, theta, Rd, level2$PPFD)
	head(level2$Af)

	#timestep change in PPFD
	level2$dPPFD <- c(NA, level2$PPFD[2:nrow(level2)] - level2$PPFD[1:(nrow(level2) - 1)])

	#identify if stepping up or not
	level2$PPFD.TF <- level2$PPFD > 0
	level2$dPPFD.TF <- level2$dPPFD != 0

	#duration of timestep and rules for carrying out operations
	level2$t.since.dPPFD <- rep(NA, nrow(level2))
	for (i in 2:nrow(level2)){
		level2$t.since.dPPFD[i] <- ifelse(level2$PPFD.TF[i],
																			ifelse(level2$PPFD.TF[i] & level2$dPPFD.TF[i], 57.6, level2$dPPFD.TF[i - 1] + 57.6),
																			0
																			)
		}

	level2$t.since.dPPFD.TF <- level2$t.since.dPPFD != 0 & !is.na(level2$t.since.dPPFD)
	level2$t.since.dPPFD.notmax <- c(NA,
																		level2$t.since.dPPFD.TF[c(1:(nrow(level2) - 1))] < level2$t.since.dPPFD.TF[c(2:nrow(level2))]
																		)
	level2$t.since.dPPFD.mid <- level2$t.since.dPPFD.notmax & level2$t.since.dPPFD.TF
	level2$dur.interval.s <- ifelse(level2$t.since.dPPFD.mid, 0, level2$t.since.dPPFD)
	level2$do.calc <- c("TRUE", level2[c(2:nrow(level2)), "dur.interval.s"] > 0)

	#initial state for start of each timestep
	#both up and down possibilities are predicted from previous timestep, then correct direction (up or down) selected based on Af versus Ai for the preceeding step
	level2$Ai.up <- c(rep(0, nrow(level2)))
	level2$Ai.down <- c(rep(0, nrow(level2)))
	level2$Ai <- c(rep(0, nrow(level2)))
	
	for (i in 2:nrow(level2)){

		L2.Afp = level2[i-1, "Af"]
		L2.t = i - (level2[i, "t.since.dPPFD"] / 57.6)

		level2[i, "Ai.up"] <- L2.Afp - (L2.Afp - level2[L2.t, "Ai"]) * exp(-level2[L2.t, "t.since.dPPFD"] / tau.up)

		level2[i, "Ai.down"] <- L2.Afp - (L2.Afp - level2[L2.t, "Ai"]) * exp(-level2[L2.t, "t.since.dPPFD"] / tau.down)

		level2[i, "Ai"] <- ifelse(level2[L2.t, "Af"] > level2[L2.t, "Ai"],
															level2[i, "Ai.up"],
															level2[i, "Ai.down"]
															)

	 }

	#establish which timesteps show a difference between instantaneous and limited kinetics
	level2$Af.Ai <- rep(0, nrow(level2))

	for (i in 1:nrow(level2)){
		if (level2[i, "do.calc"] == TRUE){
			level2[i, "Af.Ai"] <- level2[i, "Af"] - level2[i, "Ai"]
		 }
	 }

	#integrated gross assimilation

	#using rubisco limitation: applies only when PPFD is increasing
	level2$Ahat <- level2$Af * level2$t.since.dPPFD - level2$Af.Ai * tau.up + level2$Af.Ai * tau.up * exp(-level2$t.since.dPPFD / tau.up)

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