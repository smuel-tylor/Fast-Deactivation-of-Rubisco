#diurnal carbon assimilation models for planned Nature Plants article

rm(list=ls())
library(nlme)

setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rdata")
load("082005NaturePlantsnlmeLightResponses.Rdata")
load("082005NaturePlantsnlmeVcmax.Rdata")
load("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_Biochemistry/Cowpea Time-Lapse Assays/Rmodelfitting/082005NaturePlantsnlmeActivationState.Rdata")

#source code for diurnal simulation
source("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rscripts/082005predictdiurnal.R")

#need to run models with alternative parameterisations for the four genotypes
#1) full Vcmax induction parameterisation
#2) full AS parameterisation
#3) Vcmax down, AS up
#4) AS down, Vcmax up

#start by bringing together all the alternatives from the different models
objects()

AQ.fixed
AS.fixed
Vcind.fixed

d_in <- data.frame(	geno = c("Vadenantha",  "TVNu-1948",  "IT82E-16",  "IT86D-1010"),
										phi = AQ.fixed[grep("phi", rownames(AQ.fixed)), "Est"],
										Asat = AQ.fixed[grep("Asat", rownames(AQ.fixed)), "Est"],
										theta = AQ.fixed[grep("theta", rownames(AQ.fixed)), "Est"],
										Rd = AQ.fixed[grep("Rd", rownames(AQ.fixed)), "Est"],
										tau.d.AS = AS.fixed[grep("tau.d", rownames(AS.fixed)), "Est"],
										tau.a.AS = rep(AS.fixed[grep("tau.a", rownames(AS.fixed)), "Est"], 4),
										tau.d.Vc = -1200 / log(Vcind.fixed[grep("Vcmax.i", rownames(Vcind.fixed)), "Est"] / Vcind.fixed[grep("Vcmax.f", rownames(Vcind.fixed)), "Est"]),
										tau.a.Vc = rep(Vcind.fixed[grep("tau.a", rownames(Vcind.fixed)), "Est"], 4)
										)

#light regime
level2 <- read.csv("level2.csv")
#limit this to daylight hours because there's a tail of small values after sunset
#this does include 2 intervals with 0 PPFD at the start of the day that are needed to stop the algroiithm tripping up
level2 <- level2[c(283:1214),]


do.diurnal <- function(pars,light,down,up){
	
	dmods  <- lapply(pars$geno,identity)
	names(dmods) <- pars$geno
	
	for (i in 1:nrow(pars)){
		
		print(as.character(pars[i, "geno"]))
		pari = rep(NA, 6)
		names(pari) <- c("phi", "Asat", "theta", "Rd", "tau.up", "tau.down")
		pari[c("phi", "Asat", "theta", "Rd")] <- pars[i,c("phi", "Asat", "theta", "Rd")]
		pari["tau.up"]<-pars[i,up]
		pari["tau.down"]<-pars[i,down]
		pari <- unlist(pari)
		#print(pari)
		dmods[[i]] <- data.frame(geno = rep(pars[i,"geno"], nrow(level2)),
																diurnal.mod(light, pari)[[1]]
																)
		}
		
	dmods
	
}

summary.dmods <- function(dmod){
		cum.nolag <- sum(dmod[ , "Aft"], na.rm = TRUE) / 1000 #mmol/m2d 
		cum.Rubisco <- sum(dmod[ , "predict"], na.rm = TRUE) / 1000 #mmol/m2d
		perc.diff <- 100 * ((cum.nolag/cum.Rubisco) - 1)
		
		c(cum.nolag = cum.nolag,
			cum.Rubisco = cum.Rubisco,
			foregone = cum.nolag - cum.Rubisco,
			perc.diff = perc.diff
			)
	}

dmods.Vc  <- do.diurnal(d_in, level2, "tau.d.Vc", "tau.a.Vc")
summ.dmods.Vc <- sapply(dmods.Vc, summary.dmods)

dmods.AS  <- do.diurnal(d_in, level2, "tau.d.AS", "tau.a.AS")
summ.dmods.AS <- sapply(dmods.AS, summary.dmods)

dmods.VA  <- do.diurnal(d_in, level2, "tau.d.Vc", "tau.a.AS")
summ.dmods.VA <- sapply(dmods.VA, summary.dmods)

dmods.AV  <- do.diurnal(d_in, level2, "tau.d.AS", "tau.a.Vc")
summ.dmods.AV <- sapply(dmods.AV, summary.dmods)

summ.dmods <- data.frame(model = rep(c("Vc", "AS", "VA", "AV"), each = nrow(d_in)),
													geno = rep(d_in$geno, nrow(d_in)),
													t(cbind(summ.dmods.Vc, summ.dmods.AS, summ.dmods.VA, summ.dmods.AV))
													)
summ.dmods$model <- factor(summ.dmods$model, levels = c("Vc", "AS", "VA", "AV"))

write.csv(summ.dmods,"082005NaturePlantsdiurnalmodels.csv")

#Getting means and CI across genotypes

fixefci <- function(model, nlev){
	me.ci <- apply(intervals(model)$fixed[ , c(1, 2)], 1, diff)
	me.fixef <- c(fixef(model)[1], fixef(model)[1] + fixef(model)[c(2:nlev)])
	me.preds <- cbind(me.fixef, me.fixef + me.ci[c(1:nlev)] %*% cbind(-1, 1))
	me.preds <- as.data.frame(me.preds)		
	names(me.preds) <- c("Est", "lower", "upper")
	list(predictions = me.preds, CI95 = me.ci)
	}

me.diurnal <- lme(perc.diff ~ model, random = ~1 | geno, data = summ.dmods)
fixefci(me.diurnal, nlev = 4)

me.diurnal2 <- lme(foregone ~ model, random = ~1 | geno, data = summ.dmods)
fixefci(me.diurnal2, nlev = 4)

#to compare cumulatives, need to treat cum.nolag as a model for comparison
cum1 <- summ.dmods[c(1:4), c(1:3)]
cum2 <- summ.dmods[ , c(1,2,4)]
names(cum1)[3] <- names(cum2)[3] <- "cum"
cum1$model <- "absnolag"
cum<-rbind(cum1, cum2)

me.diurnal3 <- lme(cum ~ model, random = ~1 | geno, data = cum)
fixefci(me.diurnal3, nlev = 5)

save.image("082005NaturePlantsdiurnalmodels.Rdata")
