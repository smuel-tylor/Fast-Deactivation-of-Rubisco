#Diurnal carbon assimilation models for submitted Nature Plants article.
#0121, updated to work with R 4.x, here() and for style
# also reflects corrections made in nlmeVcmax.R, i.e. simpler model for Vcmax

library(here)
library(nlme)

load(here("output/nlmeLightResponses.Rdata"))
load(here("output/nlmeVcmax.Rdata"))
load(here("output/nlmeActivationState.Rdata"))

#source code for diurnal simulation
source(here("4_source_scripts_AQ_Aci_Vcmaxt_diurnal/predictdiurnal.R"))

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

#added 0121 and adapted 0321
Af <- Vcind.fixed[grep("Vcmax.f", rownames(Vcind.fixed)), "Est"]
Ai <- Vcind.fixed[grep("Vcmax.i", rownames(Vcind.fixed)), "Est"]
tau.d.Vc <- tau.down(
  Af = Af,
  Ai = Ai,
  shade.s = rep(1200, length(Af))
)

d_in <- data.frame(
  geno = levels(Vutl.noAdb$geno),
  phi = AQ.fixed[grep("phi", rownames(AQ.fixed)), "Est"],
  Asat = AQ.fixed[grep("Asat", rownames(AQ.fixed)), "Est"],
  theta = AQ.fixed[grep("theta", rownames(AQ.fixed)), "Est"],
  Rd = AQ.fixed[grep("Rd", rownames(AQ.fixed)), "Est"],
  tau.d.AS = AS.fixed[grep("tau.d", rownames(AS.fixed)), "Est"],
  tau.a.AS = rep(AS.fixed[grep("tau.a", rownames(AS.fixed)), "Est"], 4),
  tau.d.Vc = tau.d.Vc,
  tau.a.Vc = rep(Vcind.fixed[grep("tau.a", rownames(Vcind.fixed)), "Est"], 4)
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

dmods.Vc  <- do.diurnal(d_in, level2, "tau.d.Vc", "tau.a.Vc")
summ.dmods.Vc <- sapply(dmods.Vc, summary.dmods)

dmods.AS  <- do.diurnal(d_in, level2, "tau.d.AS", "tau.a.AS")
summ.dmods.AS <- sapply(dmods.AS, summary.dmods)

dmods.VA  <- do.diurnal(d_in, level2, "tau.d.Vc", "tau.a.AS")
summ.dmods.VA <- sapply(dmods.VA, summary.dmods)

dmods.AV  <- do.diurnal(d_in, level2, "tau.d.AS", "tau.a.Vc")
summ.dmods.AV <- sapply(dmods.AV, summary.dmods)

mnms <- c("Vc", "AS", "VA", "AV")

summ.dmods <- data.frame(
  model = rep(mnms, each = nrow(d_in)),
  geno = rep(d_in$geno, nrow(d_in)),
  t(cbind(summ.dmods.Vc,
          summ.dmods.AS,
          summ.dmods.VA,
          summ.dmods.AV
          )
    )
)
summ.dmods$model <- factor(summ.dmods$model, levels = mnms)

write.csv(summ.dmods,here("output/diurnalmodels.csv"),
          row.names = FALSE
          )

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

me.diurnal <- lme(perc.diff ~ model, random = ~1 | geno, data = summ.dmods)
fixefci(me.diurnal, nlev = 4)
anova(me.diurnal)

me.diurnal2 <- lme(foregone ~ model, random = ~1 | geno, data = summ.dmods)
fixefci(me.diurnal2, nlev = 4)
anova(me.diurnal2)

#to compare cumulatives, need to treat cum.nolag as a model for comparison
cum1 <- summ.dmods[c(1:4), c(1:3)]
cum2 <- summ.dmods[ , c(1, 2, 4)]
names(cum1)[3] <- names(cum2)[3] <- "cum"
cum1$model <- "absnolag"
cum<-rbind(cum1, cum2)

me.diurnal3 <- lme(cum ~ model, random = ~ 1 | geno, data = cum)
fixefci(me.diurnal3, nlev = 5)
anova(me.diurnal3)

save.image(here("output/diurnalmodels.Rdata"))
