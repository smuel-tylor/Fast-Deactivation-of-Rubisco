#Non-linear-mixed-effects models, obtaining coefficients needed in diurnal model
#This script: half times for ActivationState
#This is a separate dataset from the gas exchange, so factor order and naming
# is adjusted to be parallel with those analyses.
#0121 updated to allow for use with R 4.O and here()
#0121 revised the input .csv to remove the need for significant formatting work
# previously done in this script

library(here)
library(lattice)
library(nlme)

################################################################################
################################################################################
#bring in data and tidy

#setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_Biochemistry/Cowpea Time-Lapse Assays/Rmodelfitting")

Vutl <- read.csv(here("data/RIPE_20210121_cowpeaAS.csv"))
str(Vutl)
#update naming to be consistent with other analyses
#not needed because of revised dataset formatting
#Vutl$Genotype <- factor(replace(as.character(Vutl$Genotype),
#																Vutl$Genotype == "Adenantha",
#																"Vadenantha"
#																)
#												)

#set up factors, as needed in R 4.x
#make factors
g.levs <- c("V. adenantha", "TVNu-1948", "IT86D-1010", "IT82E-16")

Vutl$geno <- factor(Vutl$geno, levels = g.levs)
Vutl$block <- factor(Vutl$block)

p.levs <- c(sapply(g.levs, grep, x = unique(Vutl$plant), value = TRUE))
p.levs <- c(p.levs, recursive = TRUE)

Vutl$plant <- factor(Vutl$plant, levels = p.levs)

#code the individual curves to be uniquely identifiable
#0121 not needed because of data reformat
#Vutl$Replicate.ST <- factor(paste(Vutl$Genotype, Vutl$Replicate, sep = "_"))
#Vutl$time.s <- Vutl$time * 60
Vutl$time.s <- Vutl$time * 60

#pad out the missing values in each set of times*plant
#a vector that gives what would be expected in a square dataset
Ts <- unique(Vutl$time.s)
pl <- unique(Vutl$plant)
E.Tspl <- paste(rep(pl, each = length(Ts)),
                rep(Ts, times = length(pl)),
                sep = "_"
                )

Vutl2<-data.frame(
  geno = factor(
    sapply(strsplit(E.Tspl, "_"),
           function(.){ .[1] }
    ),
    levels = g.levs
  ),
  plant = factor(
    sapply(strsplit(E.Tspl, "_"),
           function(.){ paste(.[c(1, 2)], collapse = "_") }
    ),
    levels = p.levs
  ),
  time.s = as.numeric(
    sapply(strsplit(E.Tspl, "_"),
           function(.){ .[3] }
    )
  )
)

nrow(Vutl)
nrow(Vutl2)
nrow(merge(Vutl2, Vutl, all.x = TRUE))
#there are extra rows...
#must be some samples where my replicate designation covers two or more data

table(merge(Vutl2, Vutl,all.x = TRUE)$plant)
table(Vutl2$plant)
table(Vutl$plant)

#remove duplicate cases
#This was cross checked with Emmanuel Gonzalez-Escobar

#V. adenantha data
Vutl[Vutl$plant == "V. adenantha_C", ]
table(Vutl[Vutl$plant == "V. adenantha_C", "sampleID"])
#SampleID shows there are cases where the sample was apparently extracted twice
# resulting in quite different values!
Va.samps <- table(Vutl[Vutl$geno == "V. adenantha", "sampleID"])
dups <- names(Va.samps[Va.samps >= 2])
Va.dups <- Vutl[Vutl$geno == "V. adenantha" & Vutl$sampleID %in% dups, ]
Va.dups <- Va.dups[order(Va.dups$sampleID, Va.dups$extract), ]
Va.dups <- Va.dups[seq(2, nrow(Va.dups), 2), ]
rm.Va <- row.names(Va.dups)
Vutl <- Vutl[!row.names(Vutl) %in% rm.Va, ]
table(merge(Vutl2, Vutl, all.x = TRUE)$plant)
#highlights an issue with one of the TVNu
Vutl[Vutl$plant == "TVNu-1948_C", ]
#this is a simple duplication in the original file
rownames(Vutl) <- c(1:nrow(Vutl))
Vutl[Vutl$plant == "TVNu-1948_C", ]
Vutl <- Vutl[rownames(Vutl) != 212, ]
table(merge(Vutl2, Vutl, all.x = TRUE)$plant)

Vutl <- merge(Vutl2, Vutl, all.x = TRUE)

#check relevant factors
str(Vutl)
Vutl$geno
# look OK, and dataframe still ordered by geno

#looksee
xyplot(ActivationState ~ time.s | plant, data = Vutl)

##the below defunct in R 4.x,
##fix all the relevant factors by exporting this dataset and re-reading
##not actually different from previously dated versions of the same file, but changing the name here to avoid confusion
##write.csv(Vutl, "Vu_timelapse_combined_082005_square.csv", row.names = FALSE)
##Vutl <- read.csv("Vu_timelapse_combined_082005_square.csv")

##update genotype order
##Vutl$Genotype <- factor(Vutl$Genotype,
##												levels = c("Vadenantha", "TVNu-1948", "IT82E-16", "IT86D-1010")
##												)
##to make everything run smoothly, the dataframe now needs to be reordered to match this!
##Vutl<-Vutl[order(Vutl$Genotype), ]

################################################################################
################################################################################
#Analysis

head(Vutl)

#switching factors for switch-point model
Vutl$s.0 <- ifelse(Vutl$time.s < 0,  1,  0)
Vutl$s.d <- ifelse(Vutl$time.s < 0 | Vutl$time.s >= 1200,  0,  1)
Vutl$s.a <- ifelse(Vutl$time.s >= 1200,  1,  0)

#overparameterised model
Vutl.nlsList <- nlsList(
  ActivationState ~ s.0 * ASf +
    s.d *	(ASi - (ASi - ASf) * exp(-time.s / tau.d)) +
    s.a * (ASf - (ASf - ASi) * exp(-(time.s - 1200)/tau.a)) |
    plant, 
  start = c(ASi = 60, ASf = 80, tau.d = 10, tau.a = 120), 
  data = Vutl,
  na.action = na.exclude
  )

Vutl.nlsList
#one model fails

#tried different values for 'start'.
#none worked with tau.d, but this did
nls(ActivationState ~ s.0 * ASf +
											s.a * (ASf - (ASf - ASi) * exp(-(time.s - 1200) / tau.a)), 
		start = c(ASi = 60, ASf = 80, tau.a = 150), 
		data = Vutl[Vutl$plant == "V. adenantha_B", ],
		na.action = na.exclude,
		trace = TRUE
		)
#Could only be included if there was evidence
# that tau.d should not be a genotype-level parameter
Vutl.noAdb <- Vutl[Vutl$plant != "V. adenantha_B", ]

Vutl.nlsList2 <- nlsList(
  ActivationState ~ s.0 * ASf +
    s.d * (ASi - (ASi - ASf) * exp(-time.s / tau.d)) +
    s.a * (ASf - (ASf - ASi)*exp(-(time.s - 1200) / tau.a)) |
    plant, 
  start = c(ASi = 60, ASf = 80, tau.d = 10, tau.a = 120), 
  data = Vutl.noAdb,
  na.action = na.exclude
)

summary(Vutl.nlsList2)
plot(Vutl.nlsList2,  plant ~ resid(.),  abline  =  0 )
plot(intervals(Vutl.nlsList2))
#half times look difficult to identify except as an overall fixed effect
#though possibly some variation in tau.d?
#ASi and ASf look like they need modelling by genotype

#Presumably,  some level of pooling is needed
# to keep residuals homoskedastic and normal
Vutl.nls2 <- nls(
  ActivationState ~ s.0 * ASf +
    s.d * (ASi - (ASi - ASf) * exp(-time.s / tau.d)) +
    s.a * (ASf - (ASf - ASi)*exp(-(time.s - 1200) / tau.a)), 
  start = c(ASi = 60, ASf = 80, tau.d = 30, tau.a = 120), 
  data = Vutl.noAdb,
  na.action = na.exclude
)
summary(Vutl.nls2)
#residual is twice that for the fully parameterised model
plot(Vutl.nls2,  plant ~ resid(.),  abline  =  0 )
#A genotype level fixef for some params is clearly important

#Vutl.nlme2 <- nlme(
#  Vutl.nlsList2,
#  control = list(maxIter = 50,  msMaxIter  =  100,  tolerance = 5e-3)
#)
#tried a variety of settings in the above. but this won't converge.
#So, try what I expect will be a better model a priori
Vutl.nlme2 <- nlme(
  Vutl.nlsList2,
  fixed = list(ASi + ASf + tau.d + tau.a ~ 1),
  random = ASi + ASf ~ 1
  )

plot(ranef(Vutl.nlme2))
#clear evidence of genotype effects

Vutl.nlme3 <- update(
  Vutl.nlme2, 
  fixed = list(ASi + ASf + tau.d + tau.a ~ geno), 
  start = c(fixef(Vutl.nlme2)[1], 0, 0, 0,
            fixef(Vutl.nlme2)[2], 0, 0, 0,
            fixef(Vutl.nlme2)[3], 0, 0, 0,
            fixef(Vutl.nlme2)[4], 0, 0, 0), 
  random = ASi + ASf ~ 1
)

plot(ranef(Vutl.nlme3))
#genotype effects in the ranef are eliminated
anova(Vutl.nlme3)
#tau.a looks ripe for being dropped as a fixef

Vutl.nlme4 <- update(
  Vutl.nlme3, 
  fixed = list(ASi + ASf + tau.d ~ geno, tau.a ~ 1), 
  start = c(fixef(Vutl.nlme2)[1], 0, 0, 0,
            fixef(Vutl.nlme2)[2], 0, 0, 0,
            fixef(Vutl.nlme2)[3], 0, 0, 0,
            fixef(Vutl.nlme2)[4]), 
  random = ASi + ASf ~ 1
)
#supports significant differences between genotypes in initial and final values,
# as well as rate of decrease.
#But, no difference in the rate of increase in the light...
plot(ranef(Vutl.nlme4))
#still OK
anova(Vutl.nlme4)
#all significant

#these three models won't work - commented so script can be sourced

#Vutl.nlme5 <- update(
#  Vutl.nlme4, 
#  random = ASi + ASf + tau.a ~ 1
#)

#Vutl.nlme5 <- update(
#  Vutl.nlme4, 
#  random = ASi + ASf + tau.d ~ 1
#)

#Vutl.nlme5 <- update(
#  Vutl.nlme4, 
#  random = ASi + ASf + tau.a + tau.d ~ 1
#)

summary(Vutl.nlsList2)
summary(Vutl.nls2)
summary(Vutl.nlme2)
summary(Vutl.nlme3)
summary(Vutl.nlme4)
#residual for all of these is slightly larger than for the complete model
anova(Vutl.nlme2, Vutl.nlme3)
#fixed effect by genotype is clearly needed
anova(Vutl.nlme4, Vutl.nlme3)
anova(Vutl.nlme4, Vutl.nlme2)
#but tau.a genotype effect is clearly not improving the model significantly

################################################################################
################################################################################
#explore and generate outputs from the best fitting model

plot(Vutl.nlme4, plant~resid(.))
anova(Vutl.nlme4)
intervals(Vutl.nlme4)

################################################################################
#produce a nice plot of this,
# since the augPred function seems not to work for this one
t.smth <- c(-120:2400)
s.0.smth <- ifelse(t.smth < 0,  1,  0)
s.d.smth <- ifelse(t.smth < 0 | t.smth >= 1200,  0,  1)
s.a.smth <- ifelse(t.smth >= 1200,  1,  0)
p.smth <- as.character(unique(Vutl[Vutl$plant != "V. adenantha_B", "plant"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.AS <- data.frame(
  time.s = rep(t.smth, length(p.smth)), 
  s.0 = rep(s.0.smth, length(p.smth)), 
  s.d = rep(s.d.smth, length(p.smth)), 
  s.a = rep(s.a.smth, length(p.smth)), 
  geno = rep(g.smth, each = length(t.smth)), 
  plant = rep(p.smth, each = length(t.smth))
)
nd.AS$plant <- factor(nd.AS$plant, levels = levels(Vutl$plant))
nd.AS$geno <- factor(nd.AS$geno, levels = levels(Vutl$geno))

pred.Vutl.nlme4 <- data.frame(
  AS.geno = predict(Vutl.nlme4, newdata = nd.AS, level = 0), 
  AS.plant = predict(Vutl.nlme4, newdata = nd.AS, level = 1), 
  plant = nd.AS$plant,
  geno = nd.AS$geno, 
  time.s = nd.AS$time.s 
)

fhds <- c("plant", "time.s", "ActivationState")
facs.Vutl.nlme4 <- Vutl[Vutl$plant != "V. adenantha_B", fhds]
pred.Vutl.nlme4 <- merge(pred.Vutl.nlme4, facs.Vutl.nlme4, all.x = TRUE)

get_rep_list <-	function(geno, p.all){
  g.all <- p.all[p.all$geno == geno, ]
  levels(g.all$plant) <- replace(
    levels(g.all$plant),
    !levels(g.all$plant) %in% unique(g.all$plant),
    NA
  )
  by(g.all, g.all$plant, identity)
}
	
plot_rep <-	function(p.rep){
  ordp <- p.rep[order(p.rep$time.s), ]
  plot(1, 1,
       xlim = c(-120, 2400),
       ylim = c(20, 100),
       xlab = "time from start of shade (s)",
       ylab = expression(italic(S)*" (%)"),
       main = ordp$plant[1],
       type = "n",
       axes = FALSE
  )
  axis(side = 1,
       at = seq(0, 2400, 600),
       las = 1)#, labels = rep("", length(seq(0, 2400, 300))))
  axis(side = 2,
       at = seq(0, 100, 20),
       las = 1)#, labels = rep("", length(seq(0, 100, 20))))
  points(ordp$ActivationState ~ ordp$time.s)
  lines(ordp$AS.geno ~ ordp$time.s, lty = 1)
  lines(ordp$AS.plant ~ ordp$time.s, lty = 2)
  box()
}

allrep_plot <- function(geno, p.all){
  p.reps <- get_rep_list(geno, p.all)
  par(mfrow = c(2, 2))
  lapply(p.reps, plot_rep)
}

pdf(here("output/082005NaturePlantsnlmeActivationState_plot_nlme4.pdf"),
    w = 6.5, h = 6.5
    )

par(mar = c(4.5, 5.5, 3, 1), tcl = 0.4, oma = c(0, 0, 0, 0))

glist <- levels(pred.Vutl.nlme4$geno)
lapply(glist, allrep_plot, p.all = pred.Vutl.nlme4)

dev.off()

#looks OK,  and appears to justify the random effects being included
#tau.d and the activation state differ among genotypes,  tau.a is consistent

################################################################################
#fixed effects CIs as one-tailed value

cis <- apply(intervals(Vutl.nlme4)$fixed[,c(1,2)], 1, diff)
#fixed effects estimates for ASi combined with cis
fixASi <- c(fixef(Vutl.nlme4)[1],
            fixef(Vutl.nlme4)[1] + fixef(Vutl.nlme4)[c(2:4)]
            )
#checkit
cbind(fixASi, fixASi + cis[c(1:4)] %*% cbind(-1, 1))

#fixed effects estimates for ASf combined with cis
fixASf <- c(fixef(Vutl.nlme4)[5],
            fixef(Vutl.nlme4)[5] + fixef(Vutl.nlme4)[c(6:8)]
            )
#checkit
cbind(fixASf, fixASf + cis[c(5:8)] %*% cbind(-1, 1))

#fixed effects estimates for taud combined with cis
fixtaud <- c(fixef(Vutl.nlme4)[9],
             fixef(Vutl.nlme4)[9] + fixef(Vutl.nlme4)[c(10:12)]
             )
cbind(fixtaud, fixtaud + cis[c(9:12)] %*% cbind(-1, 1))

#put all of these together
AS.fixed <- rbind(
  cbind(fixASi, fixASi + cis[c(1:4)] %*% cbind(-1, 1)),
  cbind(fixASf, fixASf + cis[c(5:8)] %*% cbind(-1, 1)),
  cbind(fixtaud, fixtaud + cis[c(9:12)] %*% cbind(-1, 1)),
  cbind(fixef(Vutl.nlme4)[13],
        fixef(Vutl.nlme4)[13] + cis[c(13)] %*% cbind(-1, 1)
        )
)
AS.fixed <- as.data.frame(AS.fixed)		
names(AS.fixed)<-c("Est","lower","upper")
AS.fixed

################################################################################
#save objects

save(Vutl.nlme4,
     pred.Vutl.nlme4,
     AS.fixed,
     Vutl.noAdb,
     file=here("output/082005NaturePlantsnlmeActivationState.Rdata")
)

