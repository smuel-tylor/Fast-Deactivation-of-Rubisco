#Diurnal carbon assimilation models run as part of response to review
#for submitted Nature Plants article.
#This script specifically checks outcomes when input is nlsList outcome
#i.e., over-parameterised model run at the level of genotype means
#after summarizing those means

library(here)
library(nlme)

load(here("output/082005NaturePlantsnlmeLightResponses.Rdata"))
load(here("output/082005NaturePlantsnlmeVcmax.Rdata"))
load(here("output/082005NaturePlantsnlmeActivationState.Rdata"))

#source code for diurnal simulation
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/082005predictdiurnal.R"))

#need to run models with alternative parameterisations for the four genotypes
#1) full Vcmax induction parameterisation
#2) full AS parameterisation
#3) Vcmax down, AS up
#4) AS down, Vcmax up

#start by summarizing values from the nlsList fits
objects()

Vutl.if <- coef(Vutl.nlsList2)
Vutl.if$geno <- sapply(strsplit(row.names(Vutl.if), "_"), function(.){ .[1]})
#ensure ordering is consistent with other objects
Vutl.if$geno <- factor(Vutl.if$geno, levels = levels(Vutl.noAdb$geno))
Vutl.if.mn <- aggregate(Vutl.if[ , grep("geno", names(Vutl.if), invert = TRUE)],
                        list(geno = Vutl.if$geno),
                        mean
                        )

ind.if <- coef(ind.nlsList2)
ind.if$geno <- sapply(strsplit(row.names(ind.if), "_"), function(.){ .[1]})
#ensure ordering is consistent with other objects
ind.if$geno <- factor(ind.if$geno, levels = levels(ilGKg.df$geno))
ind.if.mn <- aggregate(ind.if[ , grep("geno", names(ind.if), invert = TRUE)],
                        list(geno = ind.if$geno),
                        mean
)

AQ.fixed
Vutl.if.mn
ind.if.mn

#calculate tau.down values for Vcmax
tau.d.Vc.if <- tau.down(
  Af = ind.if.mn$Vcmax.f,
  Ai = ind.if.mn$Vcmax.i,
  shade.s = 1200
)

d_in.if <- data.frame(
  geno = levels(Vutl.noAdb$geno),
  phi = AQ.fixed[grep("phi", rownames(AQ.fixed)), "Est"],
  Asat = AQ.fixed[grep("Asat", rownames(AQ.fixed)), "Est"],
  theta = AQ.fixed[grep("theta", rownames(AQ.fixed)), "Est"],
  Rd = AQ.fixed[grep("Rd", rownames(AQ.fixed)), "Est"],
  tau.d.AS = Vutl.if.mn$tau.d,
  tau.a.AS = Vutl.if.mn$tau.a,
  tau.d.Vc = tau.d.Vc,
  tau.a.Vc = ind.if.mn$tau.a
)

#light regime
level2 <- read.csv(here("data/level2.csv"))
#limit this to daylight hours because of a tail of small values after sunset
#The 2 intervals with 0 PPFD at the start of the day are needed to stop the
# algorithm tripping up
level2 <- level2[c(283:1214),]

do.diurnal <- function(pars, light, down, up){
  
  dmods  <- lapply(pars$geno, identity)
  names(dmods) <- pars$geno
  
  for (i in 1:nrow(pars)){
    
    print(as.character(pars[i, "geno"]))
    lrp <- c("phi", "Asat", "theta", "Rd")
    pari = rep(NA, 6)
    names(pari) <- c(lrp, "tau.up", "tau.down")
    pari[lrp] <- pars[i, lrp]
    pari["tau.up"] <- pars[i, up]
    pari["tau.down"] <- pars[i, down]
    pari <- unlist(pari)
    #print(pari)
    dmods[[i]] <- data.frame(
      geno = rep(pars[i, "geno"], nrow(level2)),
      diurnal.mod(light, pari)[[1]]
    )
  }
  
  dmods
  
}

summary.dmods <- function(dmod){
  cum.nolag <- sum(dmod[ , "Aft"], na.rm = TRUE) / 1000 #mmol/m2d 
  cum.Rubisco <- sum(dmod[ , "predict"], na.rm = TRUE) / 1000 #mmol/m2d
  #0321 modified in response to reviewer comments and
  #for consistency with Taylor & Long 2017
  perc.diff <- 100 * (1 - (cum.Rubisco / cum.nolag))
  
  c(cum.nolag = cum.nolag,
    cum.Rubisco = cum.Rubisco,
    foregone = cum.nolag - cum.Rubisco,
    perc.diff = perc.diff
  )
}

dmods.Vc.if  <- do.diurnal(d_in.if, level2, "tau.d.Vc", "tau.a.Vc")
summ.dmods.Vc.if <- sapply(dmods.Vc.if, summary.dmods)

dmods.AS.if  <- do.diurnal(d_in.if, level2, "tau.d.AS", "tau.a.AS")
summ.dmods.AS.if <- sapply(dmods.AS.if, summary.dmods)

dmods.VA.if  <- do.diurnal(d_in.if, level2, "tau.d.Vc", "tau.a.AS")
summ.dmods.VA.if <- sapply(dmods.VA.if, summary.dmods)

dmods.AV.if  <- do.diurnal(d_in.if, level2, "tau.d.AS", "tau.a.Vc")
summ.dmods.AV.if <- sapply(dmods.AV.if, summary.dmods)

mnms <- c("Vc", "AS", "VA", "AV")

summ.dmods.if <- data.frame(
  model = rep(mnms, each = nrow(d_in.if)),
  geno = rep(d_in.if$geno, nrow(d_in.if)),
  t(cbind(summ.dmods.Vc.if,
          summ.dmods.AS.if,
          summ.dmods.VA.if,
          summ.dmods.AV.if
          )
    )
)
summ.dmods.if$model <- factor(summ.dmods.if$model, levels = mnms)

write.csv(summ.dmods.if,
          here("output/082005NaturePlantsdiurnalmodels_nlsList.csv"),
          row.names = FALSE)

#Getting means and CI across genotypes

fixefci <- function(model, nlev){
  me.ci <- apply(intervals(model)$fixed[ , c(1, 2)], 1, diff)
  me.fixef <- c(
    fixef(model)[1],
    fixef(model)[1] + fixef(model)[c(2:nlev)]
  )
  me.preds <- cbind(
    me.fixef,
    me.fixef + me.ci[c(1:nlev)] %*% cbind(-1, 1)
  )
  me.preds <- as.data.frame(me.preds)		
  names(me.preds) <- c("Est", "lower", "upper")
  list(predictions = me.preds, CI95 = me.ci)
}

me.diurnal.if <- lme(perc.diff ~ model, random = ~1 | geno, data = summ.dmods.if)
fixefci(me.diurnal.if, nlev = 4)

me.diurnal.if2 <- lme(foregone ~ model, random = ~1 | geno, data = summ.dmods.if)
fixefci(me.diurnal.if2, nlev = 4)

#to compare cumulatives, need to treat cum.nolag as a model for comparison
cum1.if <- summ.dmods.if[c(1:4), c(1:3)]
cum2.if <- summ.dmods.if[ , c(1, 2, 4)]
names(cum1.if)[3] <- names(cum2.if)[3] <- "cum"
cum1.if$model <- "absnolag"
cum.if<-rbind(cum1.if, cum2.if)

me.diurnal.if3 <- lme(cum ~ model, random = ~1 | geno, data = cum.if)
fixefci(me.diurnal.if3, nlev = 5)

save.image(here("output/082005NaturePlantsdiurnalmodels_nlsList.Rdata"))
