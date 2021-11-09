#Non-linear-mixed-effects models to produce coefficients for the diurnal model
# - light response curves
# - using the same data/factors etc. as summarize.R
#updated 0121 to use R 4.x and here()

#clear workspace
#rm(list = ls())

library(here)
library(lattice)
library(nlme)

load(here("output/summarize.Rdata"))
objects()

#check factors in AQ
AQ$geno
AQ$plant

#To get the plots etc. in a nice order
AQ <- AQ[order(AQ$plant), ]

################################################################################
#Analysis

#maximal, i.e. fully parameterised model
AQ.nlsList <- nlsList(
  A ~ (
    (
      phi * Qin + Asat -
       sqrt((phi * Qin+Asat)^2 - (4 * theta * phi * Qin * Asat))
    ) /
      (2 * theta)
  ) - Rd | plant, 
  start = c(phi = 0.05, Asat = 30, theta = 0.75, Rd = 1.9),
  #Below is a quick fix for issues with some of the other columns
  data = AQ[, c("Qin", "A", "plant")],
  na.action = na.omit
)

summary(AQ.nlsList)
plot(AQ.nlsList, plant ~ resid(.), abline  =  0)
plot(intervals(AQ.nlsList))
#nothing looks very genotypey... maybe Rd?

#What is the improvement due to considering individual leaves?
#Model that ignores this:
AQ.nls <- nls(
  A ~ (
    (
      phi * Qin + Asat -
        sqrt((phi * Qin+Asat)^2 - (4 * theta * phi * Qin * Asat))
    ) / (2 * theta)
  ) - Rd, 
  start = c(phi = 0.05, Asat = 30, theta = 0.75, Rd = 1.9), 
  data = AQ,
  na.action = na.exclude
)
summary(AQ.nls)
#residual is much larger than for the fully parameterised model
plot(AQ.nls, plant ~ resid(.),  abline  =  0)
#and there is a lot of mis-fitting

AQ.nlme <- nlme(AQ.nlsList)
AQ.nlme
#Unsurprisingly, based on the nlsList object,
# random effects are reasonably independent,
# but I need to include fixed effects of genotype
AQ.nlme2 <- update(
  AQ.nlme,
  fixed = phi + Asat + theta + Rd ~ geno,
  start = c(fixef(AQ.nlme)[1], 0, 0, 0,
            fixef(AQ.nlme)[2], 0, 0, 0,
            fixef(AQ.nlme)[3], 0, 0, 0,
            fixef(AQ.nlme)[4], 0, 0, 0
  ),
  data = AQ #needed so that 'geno' is found
)
anova(AQ.nlme, AQ.nlme2)
#improves the model
anova(AQ.nlme2)
#all the fixef seem to be important, except maybe theta...
AQ.nlme2
#phi and theta could each be dropped as ranef
AQ.nlme3 <- update(
  AQ.nlme,
  fixed = list(phi + Asat + Rd ~ geno, theta ~ 1),
  start = c(fixef(AQ.nlme)[1], 0, 0, 0,
            fixef(AQ.nlme)[2], 0, 0, 0,
            fixef(AQ.nlme)[4], 0, 0, 0,
            fixef(AQ.nlme)[3]
  ),
  data = AQ #needed so that 'geno' is found
)
anova(AQ.nlme3, AQ.nlme2)
#I'm still not convinced it can be dropped...
AQ.nlme4 <- update(AQ.nlme2,
										random = Asat + theta + Rd ~ 1
										)
										
AQ.nlme5 <- update(AQ.nlme2,
										random = phi + Asat + Rd ~ 1
										)

anova(AQ.nlme2,AQ.nlme4,AQ.nlme5)
#dropping these is not useful...
#So, the model requires both fixed and random effects for all parameters

plot(AQ.nlme2)
#not ideal wrt the initial slope,
#perhaps this is not as linear at the non-rect-hyper requires
intervals(AQ.nlme2)									
#this plot failed
#plot(augPred(AQ.nlme2, primary = ~Qin, level = 0:1))

################################################################################
#produce a plant-by-plant plot matching Vcmax and ActivationState
#0121 a few changes made here in terms of specifying factors
# to make sure that objects were all consistently ordered 
q.smth <- c(0:2000)
p.smth <- as.character(unique(AQ[ , "plant"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.AQ <- data.frame(
  Qin = rep(q.smth, length(p.smth)), 
  geno = factor(
    rep(g.smth, each = length(q.smth)),
    levels = levels(AQ$geno)
  ),
  plant = factor(
    rep(p.smth, each = length(q.smth)),
    levels = levels(AQ$plant)
  )
)

pred.AQ.nlme2 <- data.frame(
  A.geno = predict(AQ.nlme2, newdata = nd.AQ, level = 0),
  A.plant = predict(AQ.nlme2, newdata = nd.AQ, level = 1),
  plant = factor(nd.AQ$plant, levels = levels(nd.AQ$plant)),
  geno = factor(nd.AQ$geno, levels = levels(nd.AQ$geno)),
  Qin = nd.AQ$Qin 
)

fhds <- c("plant", "Qin", "A")
facs.AQ.nlme2 <- AQ[, fhds]
pred.AQ.nlme2 <- merge(pred.AQ.nlme2, facs.AQ.nlme2, all.x = TRUE)

get_rep_list <-	function(geno, p.all){
  g.all <- p.all[p.all$geno == geno, ]
  levels(g.all$plant) <- replace(levels(g.all$plant),
                                 !levels(g.all$plant) %in% unique(g.all$plant),
                                 NA
  )
  by(g.all, g.all$plant, identity)
}

plot_rep <-	function(p.rep){
  ordp <- p.rep[order(p.rep$Qin), ]
  plot(1, 1,
       xlim = c(0, 2000),
       ylim = c(-3, 43),
       xlab = expression("incident PPFD " * (mu*mol~~m^-2~~s^-1)),
       ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1)),
       main = ordp$plant[1],
       type = "n",
       axes = FALSE)
  axis(side = 1, at = seq(0, 2000, 500), las = 1)
  axis(side = 2, at = seq(0, 40, 10), las = 1)
  points(ordp$A ~ ordp$Qin)
  lines(ordp$A.geno ~ ordp$Qin, lty = 1)
  lines(ordp$A.plant ~ ordp$Qin, lty = 2)
  box()
}

allrep_plot <- function(geno, p.all){
	p.reps <- get_rep_list(geno, p.all)
	par(mfrow = c(3, 2))
	lapply(p.reps, plot_rep)
	}

pdf(here("output/nlmeLightResponses_nlme2.pdf"),
    w = 6, h = 9
    )							

par(
  mar = c(4.5, 5.5, 3, 1),
  tcl = 0.4,
  oma = c(0, 0, 0, 0),
  cex = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2
)

glist <- levels(pred.AQ.nlme2$geno)
lapply(glist, allrep_plot, p.all = pred.AQ.nlme2)

dev.off()

################################################################################
#fixed effects CIs as one-tailed value

cis <- apply(intervals(AQ.nlme2)$fixed[,c(1,2)], 1, diff)

fixphi <- c(
  fixef(AQ.nlme2)[1],
  fixef(AQ.nlme2)[1] + fixef(AQ.nlme2)[c(2:4)]
)
cbind(fixphi, fixphi + cis[c(1:4)] %*% cbind(-1, 1))

fixAsat <- c(
  fixef(AQ.nlme2)[5],
  fixef(AQ.nlme2)[5] + fixef(AQ.nlme2)[c(6:8)]
)
cbind(fixAsat, fixAsat + cis[c(5:8)] %*% cbind(-1, 1))

fixtheta <- c(
  fixef(AQ.nlme2)[9],
  fixef(AQ.nlme2)[9] + fixef(AQ.nlme2)[c(10:12)]
)
cbind(fixtheta, fixtheta + cis[c(9:12)] %*% cbind(-1, 1))

fixRd <- c(
  fixef(AQ.nlme2)[13],
  fixef(AQ.nlme2)[13] + fixef(AQ.nlme2)[c(14:16)]
)
cbind(fixRd, fixRd + cis[c(13:16)] %*% cbind(-1, 1))

AQ.fixed <- rbind(
  cbind(fixphi, fixphi + cis[c(1:4)] %*% cbind(-1, 1)),
  cbind(fixAsat, fixAsat + cis[c(5:8)] %*% cbind(-1, 1)),
  cbind(fixtheta, fixtheta + cis[c(9:12)] %*% cbind(-1, 1)),
  cbind(fixRd, fixRd + cis[c(13:16)] %*% cbind(-1, 1))
)

AQ.fixed <- as.data.frame(AQ.fixed)		
names(AQ.fixed) <- c("Est", "lower", "upper")
AQ.fixed

save(AQ.fixed,
     AQ.nlme2,
     AQ,
     file = here("output/nlmeLightResponses.Rdata")
)
