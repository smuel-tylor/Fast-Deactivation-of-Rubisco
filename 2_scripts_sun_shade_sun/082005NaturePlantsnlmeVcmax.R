#Non-linear-mixed-effects models done to produce coefficients for diurnal models
#- half times for Vcmax

library(lattice)
library(nlme)
library(here)

load(here("output/082005summarize.Rdata"))
objects()

#Use "inds.list_GammastarKcogm_fixed_paired", which was identified as
# the best fitting A/ci analysis in 082005summarize.R
 
#Produce a mirror object with a simpler name
# 0121 noticed a mistake here, the object allocated was:
#  inds.list_GammastarKcoInfgm_fixed_paired.
# It should have been:
# inds.list_GammastarKcogm_fixed_paired
# This affected the nlme section and the output from that, but not the AC.j
# 0321 noticed that the lists (ilGKg and AC.fits_GammastarKcogm_fixed_paired)
#  had different orders. This needed correction in the summarize.R script
ilGKg <- inds.list_GammastarKcogm_fixed_paired
#commenting this to stop masses of console output when sourced...
#str(ilGKg)

################################################################################
################################################################################
#Identify the Vcmax limited portion of the induction
#The only effective way is A < Aj | Ci

AC.j <- function(ACisumm, ind){
	
	J = ACisumm$chosen.mod$J
	Gamma.star = ACisumm$chosen.mod$Gamma.star
	Rd = ACisumm$chosen.mod$Rd
	gm = ACisumm$chosen.mod$gm
	Pci.J = ind$Pci
	stoich = "NADPH"
					
	ind$Aj.Pci <- RuBP.limited(J, Gamma.star, Rd, gm, Pci.J, stoich)
	
	ind
}

ilGKg <- mapply(AC.j,
                AC.fits_GammastarKcogm_fixed_paired,
                ilGKg,
                SIMPLIFY = FALSE
                )

#simplify this to a dataframe, 
ilGKg.df <- data.frame(do.call(rbind, ilGKg))

#recode the timings to be zero at the start of the shade period
ilGKg.df$induction.s <- ilGKg.df$induction.s + 1200

#check/recode the columns commonly treated as factors
ilGKg.df$cult
ilGKg.df$geno
ilGKg.df$block
ilGKg.df$plant
#all look OK

#Function to compute Aj
#input here is induction data modified by AC.j to add Aj.Pci.
#This will be for one replicate
plot.AtAJ <- function(ind.AJ){
  plot(A ~ induction.s,
       data = ind.AJ,
       xlab = "time since end of shade (s)",
       ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1)),
       xlim = c(-120, 900),
       ylim = c(0, 40),
       type = "l",
       main = ind.AJ$plant[1]
       )
  lines(Aj.Pci ~ induction.s, data = ind.AJ, col = 2)
  lines(rep(60, 2), c(-1e6, 1e6), col = 5)
  lines(rep(300, 2), c(-1e6, 1e6), col = 3)
  lines(rep(600, 2), c(-1e6, 1e6), col = 4)
}

#Function looping the above.
#Here, input is AC.j modified induction data for multiple replicates
# and structured by 'geno'
plot.AtAJ.bygeno <- function(inds.AJ){
  for (i in levels(inds$geno)){
    par(mfrow = c(3, 2), mar = c(4, 5.5, 2, 0.5), las = 1, cex = 1.2)
		lapply(
		  inds.AJ[grep(i, names(inds.AJ), value = TRUE)],
		  plot.AtAJ
		  )
  }
}

#pdf for output
pdf(here("output/082005NaturePlantsnlmeVcmax_assessJlimitation.pdf"),
    paper = "a4",
    w = 8, h = 11
)

header.page(
  main = expression(
    atop("Does "*italic(A)*"(black line) intersect with",
         italic(A)[J] * " (red line),  predicted at steady state " *
           italic(J) * " & " * italic(c)[i] * " during induction"
    )
  )
)

plot.AtAJ.bygeno(ilGKg)

dev.off()

#don't know that 5 min will be enough for good fits,
#but if it works,  it's justified by these plots

#limit the timesequences to keep them consistent
# bearing in mind the 5 min post-shade window
ilGKg.df <- ilGKg.df[ilGKg.df$induction.s >= -120 &
                       ilGKg.df$induction.s <= 1200 + 300, ]

xyplot(Astar*5 + Vcmax.t ~ induction.s | plant, data = ilGKg.df)

#Retain only the columns actually needed for the model
# to ensure an na.omit below works OK (because FLR values are usually all NA)
ilGKg.df <- ilGKg.df[ , c("induction.s", "geno", "plant", "Astar", "Vcmax.t")]

#Make the dataset square so that resampling for CI will work.
#Vectors that give what would be expected in a square dataset:
Ts <- unique(ilGKg.df$induction.s)
pl <- unique(ilGKg.df$plant)
E.plTs <- paste(rep(pl, each = length(Ts)),
                rep(Ts, times = length(pl)),
                sep = "_"
                )

#make a new df using this
#revised 0121 to make sure that factors are specified parallel to ilGKg.df
ilGKg.df2 <- data.frame(
  geno = factor(
    sapply(strsplit(E.plTs, "_"),
           function(.){ .[1] }
    ),
    levels = levels(ilGKg.df$geno)
  ),
  plant = factor(
    sapply(strsplit(E.plTs, "_"),
           function(.){ paste(.[c(1, 2)], collapse = "_") }
    ),
    levels = levels(ilGKg.df$plant)
  ),
  induction.s = as.numeric(
    sapply(strsplit(E.plTs, "_"),
           function(.){ .[3] }
    )
  )
)

#checking this will work
head(ilGKg.df2)

nrow(ilGKg.df)
nrow(ilGKg.df2)

nrow(merge(ilGKg.df2, ilGKg.df, all.x = TRUE))#seems to work nicely
table(merge(ilGKg.df2, ilGKg.df, all.x = TRUE)$plant)#as above

#do a merge to replace ilGKg.df with the square data
ilGKg.df <- merge(ilGKg.df2, ilGKg.df, all.x = TRUE)
#clean up
rm(ilGKg.df2, Ts, pl, E.plTs)

#check
str(ilGKg.df)
#this looks correct

#to make everything run smoothly,  order ilGKg.df by geno
ilGKg.df <- ilGKg.df[order(ilGKg.df$geno), ]

#add a 'time' variable that includes NA where samples were not collected
ilGKg.df$inductionNA.s <- replace(ilGKg.df$induction.s,
                                  is.na(ilGKg.df$Vcmax.t),
                                  NA
                                  )

#quick plots
xyplot(Astar*5 + Vcmax.t ~ induction.s | plant, data = ilGKg.df)
#this matches the same plot on 124, so everything is parallel
xyplot(Astar ~ Vcmax.t | plant, data = ilGKg.df)
#Since Astar is assumed to be Vcmax limited and other parameters are constant,
#Astar and Vcmax.t are close to directly proportional
# with a curvature determined by the A/ci response...

#Modelling Vcmax.t is more obviously parallel with Rubisco activation state
# than modelling Astar.
#Astar so calculated is equivalent to Vcmax.t,
# and is only correct over the Vcmax limited range.
# True Astar will be smaller in J-limited range, i.e., after about 5 min,
# but extrapolation of Vcmax asymptote should still hold,
# irrespective of which factor is limiting.

################################################################################
################################################################################
#Timeseries analyses using nlme

#Noting that this is only valid for the initial post-shade phase:
# a piecewise model incorporating shade/pre-shade is not appropriate for Vcmax
# because we do not know that Vcmax is limiting prior to shade

#To make the piecewise model work,
#I need to add appropriate time sequences & time sequence factors to the data
#as noted above: unlike nlmeActivationState.R,  I cannot model the 0-1200 period
# of shade so code here extracts that data as ilGKg.df.ind
ilGKg.df$s.a <- ifelse(ilGKg.df$induction.s > 1260,  1,  0) #fails if >=1260!
ilGKg.df.ind <- ilGKg.df[ilGKg.df$s.a == 1, ]

#first model: over-parameterised, i.e., complete rep-by-rep
ind.nlsList <- nlsList(
  Vcmax.t ~ (Vcmax.f - (Vcmax.f - Vcmax.i) *
               exp( -(induction.s - 1200) / tau.a)) | plant,
  start = c(Vcmax.i = 60, Vcmax.f = 170, tau.a = 120),
  data = ilGKg.df.ind,
  na.action = na.exclude
  )

summary(ind.nlsList)
plot(ind.nlsList,  plant ~ resid(.),  abline  =  0 )
plot(intervals(ind.nlsList))
#As you might anticipate from noisy data,
# some really poorly constrained estimates of the coefs

#0321 - remove IT86D-1010_9 and V.adenantha_8 because of high residual variance
# -problematic in terms of homoskedasticity
# -also in terms of the basic model because some early data exceed the asymptote
ilGKg.df.sub <- ilGKg.df.ind[ilGKg.df.ind$plant != "IT86D-1010_9" &
                                     ilGKg.df.ind$plant != "V. adenantha_8", ]

ind.nlsList2 <- nlsList(
  Vcmax.t ~ (Vcmax.f - (Vcmax.f - Vcmax.i) *
               exp( -(induction.s - 1200) / tau.a)) | plant,
  start = c(Vcmax.i = 60, Vcmax.f = 170, tau.a = 120),
  data = ilGKg.df.sub,
  na.action = na.exclude
)

summary(ind.nlsList2)
plot(ind.nlsList2,  plant ~ resid(.),  abline  =  0 )
plot(intervals(ind.nlsList2))
#note that parameters are still not always clearly identifiable
# Vcmax.f anf tau.a are poorly constrained where induction was shallow

#Some level of pooling will help keep residuals homoskedastic and normal.
#This model is totally pooled, i.e., under-parameterised.
ind.nls <- nls(
  Vcmax.t ~ (Vcmax.f - (Vcmax.f - Vcmax.i) * exp(-(induction.s - 1200) / tau.a)),
  start = c(Vcmax.i = 60, Vcmax.f = 170, tau.a = 120),
  data = ilGKg.df.sub, na.action = na.exclude
  )

summary(ind.nls)
#As you would anticipate,
# the residual is much larger than for the fully parameterised model,
# supporting a need for random effects
plot(ind.nls,  plant ~ resid(.),  abline  =  0 )
#Explains why:
# some plants are deviating strongly from the mean 

#build a random effects model starting from the over-parameterised ind.nlsList2
#have to boost the iterations a lot here to force this to work
ind.nlme <- nlme(ind.nlsList2, control = list(maxIter = 1000, msMaxIter = 1000))
ind.nlme
#this looks OK...

plot(ranef(ind.nlme))
#No obvious structure, based on this plot
# except... maybe Vcmax.i differs between genotypes?? 
# incorporate and test for necessary fixed effects of genotype
ind.nlme2 <- update(ind.nlme,
                    fixed = Vcmax.i + Vcmax.f + tau.a ~ geno,
                    start = c(fixef(ind.nlme)[1], 0, 0, 0,
                              fixef(ind.nlme)[2], 0, 0, 0,
                              fixef(ind.nlme)[3], 0, 0, 0
                              )
                    )
#works, with a warning about non-convergence
anova(ind.nlme, ind.nlme2)
#this is marginally a bad step... AIC increases, but we need to know about genotype effects
summary(ind.nlme2)
#nothing looks useful here
anova(ind.nlme2)
#This indicates that dropping tau.a as a genotype-level fixed effect makes sense
ind.nlme3 <- update(ind.nlme,
										fixed = list(Vcmax.i + Vcmax.f ~ geno, tau.a ~ 1),
										start = c(fixef(ind.nlme)[1], 0, 0, 0,
															fixef(ind.nlme)[2], 0, 0, 0,
															fixef(ind.nlme)[3]
															)
										)
anova(ind.nlme, ind.nlme2, ind.nlme3)
#smaller than ind.nlme2, but also increasing AIC compared with ind.nlme
#Similar effect on likelihood to ind.nlme2
pairs(ind.nlme3)
anova(ind.nlme3)
#looks like all fixef are reasonable here
summary(ind.nlme3)
#but only Vcmax.i coefs look to be significant

#so check Vcmax.f.geno
ind.nlme4 <- update(ind.nlme,
                    fixed = list(Vcmax.i ~ geno, Vcmax.f + tau.a ~ 1),
                    start = c(fixef(ind.nlme)[1], 0, 0, 0,
                              fixef(ind.nlme)[2],
                              fixef(ind.nlme)[3]
                    )
)
anova(ind.nlme, ind.nlme2, ind.nlme3, ind.nlme4)
#no issue with dropping this despite low Wald test p... (F value was marginal).
# This model with only Vcmax.i, which was identified graphically to begin,
# is now better than the complete random effect model by AIC
pairs(ind.nlme4)
anova(ind.nlme4)
summary(ind.nlme4)

#belt & braces sanity check if you drop Vcmax.i instead or retain tau.a
ind.nlme5 <- update(ind.nlme,
                    fixed = list(Vcmax.f ~ geno, Vcmax.i + tau.a ~ 1),
                    start = c(fixef(ind.nlme)[2], 0, 0, 0,
                              fixef(ind.nlme)[1],
                              fixef(ind.nlme)[3]
                    )
)

anova(ind.nlme, ind.nlme2, ind.nlme3, ind.nlme4, ind.nlme5)
#ind.nlme4, is clearly the best fixef model

#Is it reasonable to drop any of the random effects?

#drop tau.a as a rep-by-rep random effect
ind.nlme6 <- update(ind.nlme4,
                    random = Vcmax.i + Vcmax.f ~ 1
)
#alternatively, drop Vcmax.f as a rep-by-rep random effect
ind.nlme7 <- update(ind.nlme4, 
                    random = Vcmax.i + tau.a ~ 1
)
#alternatively, drop Vcmax.i as a rep-by-rep random effect
ind.nlme8 <- update(ind.nlme4, 
                    random = Vcmax.f + tau.a ~ 1
)

#compare these with the model that has a full complement of random effects
# and minimised fixed effects
anova(ind.nlme4, ind.nlme6, ind.nlme7, ind.nlme8)
#None of these can be dropped as random effects
# without significant decrease in logLik, and increase in AIC

#0321 The below is modified to reflect use of ilGKg.df.sub & ind.nlme4
################################################################################
################################################################################
#Plots and output

################################################################################
#produce a plot equivalent to the augPred function
t.smth <- c(1200:1800)
p.smth <- as.character(unique(ilGKg.df.sub[ , "plant"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.Vc <- data.frame(induction.s = rep(t.smth, length(p.smth)), 
                    geno = rep(g.smth, each = length(t.smth)), 
                    plant = rep(p.smth, each = length(t.smth))
)
nd.Vc$plant <- factor(nd.Vc$plant, levels = levels(ilGKg.df.sub$plant))
nd.Vc$geno <- factor(nd.Vc$geno, levels = levels(ilGKg.df.sub$geno))

pred.ind.nlme4 <- data.frame(
  Vc.geno = predict(ind.nlme4, newdata = nd.Vc, level = 0), 
	Vc.plant = predict(ind.nlme4, newdata = nd.Vc, level = 1), 
	plant = nd.Vc$plant,
	geno = nd.Vc$geno, 
	induction.s = nd.Vc$induction.s
  )

fhds <- c("plant", "induction.s", "Vcmax.t")
facs.ind.nlme4 <- ilGKg.df.sub[, fhds]
pred.ind.nlme4 <- merge(pred.ind.nlme4, facs.ind.nlme4, all.x = TRUE)

get_rep_list <-	function(geno, p.all){
  g.all <- p.all[p.all$geno == geno, ]
  levels(g.all$plant) <- replace(levels(g.all$plant),
                                 !levels(g.all$plant) %in% unique(g.all$plant),
                                 NA)
  by(g.all, g.all$plant, identity)
}
	
plot_rep <-	function(p.rep){
  ordp <- p.rep[order(p.rep$induction.s), ]
  plot(1, 1,
       xlim = c(1200, 1800),
       ylim = c(0, 400),
       xlab = "time from end of shade (s)",
       ylab = expression(italic(V)["c,max"]~~(mu * mol~~m^-2~~s^-1)),
       main = ordp$plant[1],
       type = "n",
       axes = FALSE)
  axis(side = 1,
       at = seq(1200, 1800, 120),
       labels = seq(1200, 1800, 120) - 1200,
       las = 1
       )
  axis(side = 2,
       at = seq(0, 400, 100),
       las = 1
       )
  points(ordp$Vcmax.t ~ ordp$induction.s)
  lines(ordp$Vc.geno ~ ordp$induction.s, lty = 1)
  lines(ordp$Vc.plant ~ ordp$induction.s, lty = 2)
  box()
}

allrep_plot <- function(geno, p.all){
  p.reps <- get_rep_list(geno, p.all)
  par(mfrow = c(3, 2))
  lapply(p.reps, plot_rep)
}

pdf(here("output/082005NaturePlantsnlmeVcmax_plot_nlme4.pdf"),
    w = 6, h = 9
    )							

par(mar = c(4.5, 5.5, 3, 1),
    tcl = 0.4,
    oma = c(0, 0, 0, 0),
    cex = 1.2,
    cex.lab = 1.2,
    cex.axis = 1.2
    )

glist <- levels(pred.ind.nlme4$geno)
lapply(glist, allrep_plot, p.all = pred.ind.nlme4)

dev.off()

################################################################################
#Here's a plot of the initial part of the Vcmax ~ t response,
#showing that after the initial 5 min period when Vcmax is always limiting
# there are only very marginal effects
#0121, removed genotype effect, as this is not present in ind.nlme7
# used ind.nlme7, and renamed objects to reflect 7 instead of 6
pdf(here("output/082505NaturePlantsnlmeVcmax_extrapolation.pdf"),
    paper = "a4",
    w = 8, h = 8
    )

nd <- data.frame(induction.s = rep(seq(1200, 2400, 10), 4),
                 geno = rep(levels(ilGKg.df.sub$geno),
                            each = length(seq(1200, 2400, 10))
                 )
)
in4 <- predict(ind.nlme4, 
								newdata = nd, 
								level = 0
								)
in4 <- data.frame(nd, Vcmax.t = in4)
par(mfrow = c(1, 1), mar = c(5, 5, 5, 5), las = 1)
plot(Vcmax.t ~ induction.s,
     data = ilGKg.df.sub[ilGKg.df.sub$induction.s > 1260, ], 
     las = 1,
     pch = 21,
     col = NULL,
     bg = rgb(0, 0, 0, alpha = 0.2), 
     xlab = "time since shade (min)",
     xlim = c(1200, 1800),
     xaxt = "n", 
     ylab = expression(italic(V)["c, max"]~~(mu*mol~~m^-2~~s^-1)),
     ylim = c(0, 340),
     yaxt = "n"
)
axis(side = 1,
     at = c(1200, 1500, 1800),
     labels = c(0, 5, 10)
)
axis(side = 2,
     at = seq(0, 300, 100)
)
title(
  main = expression(
    atop(italic(V)["c, max"] * " timeseries from gas exchange in four " *
           italic(Vigna) * " genotypes",
         "effect of projection beyond " * italic(V)["c.max"] * " limited range"
         )
    )
  )

by(in4,
   in4$geno,
   function(.){
     lines(Vcmax.t ~ induction.s, data = .)
     points(30 * diff(.$Vcmax.t) ~ .$induction.s[c(2:nrow(.))],
            col = rgb(1, 0, 0, alpha = 0.3),
            pch = 21
     )
     }
)

#points(30 * diff(177 - (177 - 88) * exp(-seq(0, 600, 10) / 129)) ~
#         seq(1210, 1800, 10),
#       col = 2,
#       pch = 19
#       )

axis(side = 4,
     at = c(0, 60, 120, 180, 240, 300),
     labels = c(0, 60, 120, 180, 240, 300) / 30,
     col = 2,
     col.axis = 2,
     las = 1
)

mtext(side = 4,
      line = 3,
      expression("change in " * italic(V)["c, max"] * " per 10 s interval"),
      col = 2,
      las = 3
)

dev.off()

#Vcmax fixed effects CIs
cis <- apply(intervals(ind.nlme4, which = "fixed")$fixed[ , c(1, 2)], 1, diff)
#fixed effects estimates for ASi combined with cis
nmsfe <- names(fixef(ind.nlme4))
all(nmsfe == names(cis))

fixVci <- c(fixef(ind.nlme4)[grep("Intercept", nmsfe)],
            fixef(ind.nlme4)[grep("Intercept", nmsfe)] +
              fixef(ind.nlme4)[grep("Vcmax.i.geno", nmsfe)]
)

fixVcf <- fixef(ind.nlme4)[grep("Vcmax.f", nmsfe)]

fixtau <- fixef(ind.nlme4)[grep("tau", nmsfe)]

Vcind.fixed <- rbind(
  cbind(fixVci, fixVci + cis[grep("Vcmax.i", nmsfe)] %*% cbind(-1, 1)),
  cbind(fixVcf, fixVcf + cis[grep("Vcmax.f", nmsfe)] %*% cbind(-1, 1)),
  cbind(fixtau, fixtau + cis[grep("tau", nmsfe)] %*% cbind(-1, 1))
)
Vcind.fixed <- as.data.frame(Vcind.fixed)		
names(Vcind.fixed) <- c("Est", "lower", "upper")
Vcind.fixed

################################################################################
################################################################################
#Save outputs

#updated 0321
save(in4,
     ind.nlsList2,
     ind.nlme4,
     Vcind.fixed,
     ilGKg.df.sub,
     file = here("output/082005NaturePlantsnlmeVcmax.Rdata")
     )
