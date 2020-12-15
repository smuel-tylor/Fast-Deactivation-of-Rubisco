#Non-linear-mixed-effects model building for the coefficients needed in diurnal model
#half times for ActivationState

#The order of the "Genotype" factor is adjusted to be: Adenantha,  TVNu-1948,  IT82E-16,  IT86D-1010
#and square datasets produced for subsequent use in bootstrapping CIs

#This modifies the timesequencemodel approach to characterising parameters needed for gas exchange to work with AQ fitting
#intended to provide input for diurnal models that is consistent with
#timesequencemodel.R
#timesequencemodel_gasexchange.R

#as per the biochemistry dataset, the "Genotype" factor is adjusted at the start to be: Adenantha, TVNu-1948, IT82E-16, IT86D-1010
#the code is prone to getting mixed up around this when editing, so requires checking on completion

#################################################
#THIS CODE FAILS IN R 4.0 BECAUSE read.csv differs
##################################################

#clear workspace
rm(list=ls())

#library(readxl)
library(lattice)
library(nlme)

#to make the piecewise model work, I need to add appropriate time sequence factors to the dataset
#to be able to compute confidence intervals using bootstraps, I need the dataset to be square

#bring in data
setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_Biochemistry/Cowpea Time-Lapse Assays/Rmodelfitting")

Vutl <- read.csv("Vu_timelapse_combined_032010.csv")
str(Vutl)
#update naming to be consistent with other analyses
Vutl$Genotype <- factor(replace(as.character(Vutl$Genotype),
																Vutl$Genotype == "Adenantha",
																"Vadenantha"
																)
												)

#code the individual curves to be uniquely identifiable
Vutl$Replicate.ST <- factor(paste(Vutl$Genotype, Vutl$Replicate, sep = "_"))
Vutl$TimePoint.s <- Vutl$TimePoint * 60

#pad out the missing values in each set of TimePoints*Replicate.ST
#Vutl$Replicate.ST.s<-paste(Vutl$Replicate.ST,Vutl$TimePoint.s,sep="_")
#a vector that gives what would be expected in a square dataset
Ts <- unique(Vutl$TimePoint.s)
RST <- unique(Vutl$Replicate.ST)
E.RSTs <- paste(rep(RST, each = length(Ts)), rep(Ts, times = length(RST)), sep = "_")

Vutl2<-data.frame(Genotype = sapply(strsplit(E.RSTs, "_"),
																		function(.){ .[1] }
																		),
									Replicate.ST = sapply(strsplit(E.RSTs, "_"),
																				function(.){ paste(.[c(1, 2)], collapse = "_") }
																				),
									TimePoint.s = sapply(strsplit(E.RSTs, "_"),
																				function(.){ .[3] }
																				)
									)

nrow(Vutl)
nrow(Vutl2)
nrow(merge(Vutl2, Vutl,all.x = TRUE))#there are extra rows... must be some samples where my replicate designation covers two or more data

table(merge(Vutl2, Vutl,all.x = TRUE)$Replicate.ST)
table(Vutl2$Replicate.ST)
table(Vutl$Replicate.ST)

#remove duplicate cases

#Vadenantha data
Vutl[Vutl$Replicate.ST == "Vadenantha_C", ]
#SampleID shows these are cases where the sample was apparently extracted twice - often resulting in quite different values!
Va.samps <- table(Vutl[Vutl$Genotype == "Vadenantha", c("SampleID")])
dups <- names(Va.samps[Va.samps >= 2])
Va.dups <- Vutl[Vutl$Genotype == "Vadenantha" & Vutl$SampleID %in% dups, ]
Va.dups <- Va.dups[order(Va.dups$SampleID, Va.dups$Extract), ]
Va.dups <- Va.dups[seq(2, nrow(Va.dups), 2), ]
rm.Va <- row.names(Va.dups)
Vutl <- Vutl[!row.names(Vutl) %in% rm.Va, ]
table(merge(Vutl2, Vutl, all.x = TRUE)$Replicate.ST)

#for checking with EGE
Va.dups[rm.Va, "Extract"]

#also an issue with one of the TVNu
Vutl[Vutl$Replicate.ST == "TVNu-1948_C", ]
#this is a simple duplication in the original file
rownames(Vutl) <- c(1:nrow(Vutl))
Vutl[Vutl$Replicate.ST == "TVNu-1948_C", ]
Vutl <- Vutl[c(1:211, 213:nrow(Vutl)), ]
table(merge(Vutl2, Vutl, all.x = TRUE)$Replicate.ST)

Vutl <- merge(Vutl2, Vutl, all.x = TRUE)

#fix all the relevant factors by exporting this dataset and re-reading
#not actually different from previously dated versions of the same file, but changing the name here to avoid confusion
write.csv(Vutl, "Vu_timelapse_combined_082005_square.csv", row.names = FALSE)
Vutl <- read.csv("Vu_timelapse_combined_082005_square.csv")
str(Vutl)

#update genotype order
Vutl$Genotype <- factor(Vutl$Genotype,
												levels = c("Vadenantha", "TVNu-1948", "IT82E-16", "IT86D-1010")
												)
#to make everything run smoothly, the dataframe now needs to be reordered to match this!
Vutl<-Vutl[order(Vutl$Genotype), ]

xyplot(ActivationState ~ TimePoint.s | Replicate.ST, data = Vutl)

head(Vutl)

Vutl$s.0 <- ifelse(Vutl$TimePoint.s < 0,  1,  0)
Vutl$s.d <- ifelse(Vutl$TimePoint.s < 0 | Vutl$TimePoint.s >= 1200,  0,  1)
Vutl$s.a <- ifelse(Vutl$TimePoint.s >= 1200,  1,  0)

Vutl.nlsList <- nlsList(ActivationState ~ s.0 * ASf +
																					s.d *	(ASi - (ASi - ASf) * exp(-TimePoint.s / tau.d)) +
																					s.a * (ASf - (ASf - ASi) * exp(-(TimePoint.s - 1200)/tau.a)) |
																					Replicate.ST, 
												start = c(ASi = 60, ASf = 80, tau.d = 10, tau.a = 120), 
												data = Vutl, na.action = na.exclude)
												
xyplot(ActivationState ~ TimePoint.s | Replicate.ST, data = Vutl)

nls(ActivationState ~ s.0 * ASf +
											s.a * (ASf - (ASf - ASi) * exp(-(TimePoint.s - 1200) / tau.a)), 
		start = c(ASi = 60, ASf = 80, tau.a = 150), 
		data = Vutl[Vutl$Replicate.ST == "Vadenantha_B", ],
		na.action = na.exclude,
		trace = T
		)
#tried a variety of parameterisations,  but this Replicate won't allow tau.d to be fit,
#so can only be included if there is evidence that tau.d should be a genotype-level parameter
Vutl.noAdb <- Vutl[Vutl$Replicate.ST != "Vadenantha_B", ]

Vutl.nlsList2 <- nlsList(ActivationState ~ s.0 * ASf +
																						s.d * (ASi - (ASi - ASf) * exp(-TimePoint.s / tau.d)) +
																						s.a * (ASf - (ASf - ASi)*exp(-(TimePoint.s - 1200) / tau.a)) |
																						Replicate.ST, 
													start = c(ASi = 60, ASf = 80, tau.d = 10, tau.a = 120), 
													data = Vutl.noAdb,
													na.action = na.exclude
													)
													
summary(Vutl.nlsList2)
plot(Vutl.nlsList2,  Replicate.ST ~ resid(.),  abline  =  0 )
plot(intervals(Vutl.nlsList2))
#half times look difficult to identify except as an overall fixed effect
#though possibly some variation in tau.d?
#ASi and ASf look like they need modelling by genotype

#presumably though,  some level of pooling is needed to keep residuals homoskedastic and normal
Vutl.nls2 <- nls(ActivationState ~ s.0 * ASf +
																		s.d * (ASi - (ASi - ASf) * exp(-TimePoint.s / tau.d)) +
																		s.a * (ASf - (ASf - ASi)*exp(-(TimePoint.s - 1200) / tau.a)), 
									start = c(ASi = 60, ASf = 80, tau.d = 30, tau.a = 120), 
									data = Vutl.noAdb,
									na.action = na.exclude
									)
summary(Vutl.nls2)
#residual is twice that for the fully parameterised model
plot(Vutl.nls2,  Replicate.ST ~ resid(.),  abline  =  0 )
#explains why,  & indicates that a genotype level fixef for some params is important

#Vutl.nlme2 <- nlme(Vutl.nlsList2, control = list(maxIter = 50,  msMaxIter  =  100,  tolerance = 5e-3))
#tried a variety of settings but this won't converge
#so try what I expect will be a better model a priori
Vutl.nlme2 <- nlme(Vutl.nlsList2, 
										fixed = list(ASi + ASf + tau.d + tau.a ~ 1), 
										#start = c(Vutl.nls2[1], 0, 0, 0, Vutl.nls2[2], 0, 0, 0, Vutl.nls2[3], 0, 0, 0, Vutl.nls2[4], 0, 0, 0), 
										random = ASi + ASf ~ 1
										)
plot(ranef(Vutl.nlme2))#clear evidence of genotype effects
Vutl.nlme3 <- update(Vutl.nlme2, 
											fixed = list(ASi + ASf + tau.d + tau.a ~ Genotype), 
											start = c(fixef(Vutl.nlme2)[1], 0, 0, 0,
																fixef(Vutl.nlme2)[2], 0, 0, 0,
																fixef(Vutl.nlme2)[3], 0, 0, 0,
																fixef(Vutl.nlme2)[4], 0, 0, 0), 
											random = ASi + ASf ~ 1
											)
plot(ranef(Vutl.nlme2))#genotype effects in the ranef are gone
anova(Vutl.nlme3)
Vutl.nlme4 <- update(Vutl.nlme3, 
											fixed = list(ASi + ASf + tau.d ~ Genotype, tau.a~1), 
											start = c(fixef(Vutl.nlme2)[1], 0, 0, 0,
																fixef(Vutl.nlme2)[2], 0, 0, 0,
																fixef(Vutl.nlme2)[3], 0, 0, 0,
																fixef(Vutl.nlme2)[4]), 
											random = ASi + ASf ~ 1
											)
#supports significant differences between genotypes in initial and final values,  as well as rate of decrease
#but no difference in the rate of increase in the light...

Vutl.nlme5 <- update(Vutl.nlme4, 
											random = ASi + ASf + tau.a ~ 1
											)#won't work

Vutl.nlme5 <- update(Vutl.nlme4, 
											random = ASi + ASf + tau.d ~ 1
											)#won't work

Vutl.nlme5 <- update(Vutl.nlme4, 
											random = ASi + ASf + tau.a + tau.d ~ 1
											)#won't work

summary(Vutl.nlsList2)
summary(Vutl.nls2)
summary(Vutl.nlme2)
summary(Vutl.nlme3)
summary(Vutl.nlme4)
#residual for all of these is slightly larger than for the complete model
anova(Vutl.nlme2, Vutl.nlme3)#fixed effect by genotype is clearly needed
anova(Vutl.nlme4, Vutl.nlme3)
anova(Vutl.nlme4, Vutl.nlme2)#but tau.a genotype effect is clearly not improving the model significantly

#explore the best fitting model
plot(Vutl.nlme4, Replicate.ST~resid(.))
anova(Vutl.nlme4)
intervals(Vutl.nlme4)

#produce a nice plot of this, since the augPred function seems not to work for this one
t.smth <- c(-120:2400)
s.0.smth <- ifelse(t.smth < 0,  1,  0)
s.d.smth <- ifelse(t.smth < 0 | t.smth >= 1200,  0,  1)
s.a.smth <- ifelse(t.smth >= 1200,  1,  0)
p.smth <- as.character(unique(Vutl[Vutl$Replicate.ST != "Vadenantha_B", "Replicate.ST"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.AS <- data.frame(TimePoint.s = rep(t.smth, length(p.smth)), 
											s.0 = rep(s.0.smth, length(p.smth)), 
											s.d = rep(s.d.smth, length(p.smth)), 
											s.a = rep(s.a.smth, length(p.smth)), 
											Genotype = rep(g.smth, each = length(t.smth)), 
											Replicate.ST = rep(p.smth, each = length(t.smth))
											)
nd.AS$Replicate.ST <- factor(nd.AS$Replicate.ST, levels = levels(Vutl$Replicate.ST))
nd.AS$Genotype <- factor(nd.AS$Genotype, levels = levels(Vutl$Genotype))

pred.Vutl.nlme4 <- data.frame(AS.geno = predict(Vutl.nlme4, newdata = nd.AS, level = 0), 
															AS.plant = predict(Vutl.nlme4, newdata = nd.AS, level = 1), 
															Replicate.ST = nd.AS$Replicate.ST,
															Genotype = nd.AS$Genotype, 
															TimePoint.s = nd.AS$TimePoint.s 
															)
															
fhds <- c("Replicate.ST", "TimePoint.s", "ActivationState")
facs.Vutl.nlme4 <- Vutl[Vutl$Replicate.ST != "Vadenantha_B", fhds]
pred.Vutl.nlme4 <- merge(pred.Vutl.nlme4, facs.Vutl.nlme4, all.x = T)

get_rep_list <-	function(geno, p.all){
	g.all <- p.all[p.all$Genotype == geno, ]
	levels(g.all$Replicate.ST) <- replace(levels(g.all$Replicate.ST),
																				!levels(g.all$Replicate.ST) %in% unique(g.all$Replicate.ST),
																				NA)
	by(g.all, g.all$Replicate.ST, identity)
	}
	
plot_rep <-	function(p.rep){
	ordp <- p.rep[order(p.rep$TimePoint.s), ]
	plot(1, 1,
				xlim = c(-120, 2400),
				ylim = c(20, 100),
				xlab = "time from start of shade (s)",
				ylab = expression(italic(S)*" (%)"),
				main = ordp$Replicate.ST[1],
				type = "n",
				axes = FALSE)
	axis(side = 1, at = seq(0, 2400, 600), las = 1)#, labels = rep("", length(seq(0, 2400, 300))))
	axis(side = 2, at = seq(0, 100, 20), las = 1)#, labels = rep("", length(seq(0, 100, 20))))
	points(ordp$ActivationState ~ ordp$TimePoint.s)
	lines(ordp$AS.geno ~ ordp$TimePoint.s, lty = 1)
	lines(ordp$AS.plant ~ ordp$TimePoint.s, lty = 2)
	box()
	}

allrep_plot <- function(geno, p.all){
	p.reps <- get_rep_list(geno, p.all)
	par(mfrow = c(2, 2))
	lapply(p.reps, plot_rep)
	}

pdf("082005NaturePlantsnlmeActivationState_plot_nlme4.pdf", w = 6.5, h = 6.5)							

par(mar = c(4.5, 5.5, 3, 1), tcl = 0.4, oma = c(0, 0, 0, 0))

glist <- levels(pred.Vutl.nlme4$Genotype)
lapply(glist, allrep_plot, p.all = pred.Vutl.nlme4)

dev.off()

#looks OK,  and appears to justify the random effects being included
#tau.d and the activation state differ among genotypes,  tau.a is consistent

#fixed effects CIs as one-tailed value
cis <- apply(intervals(Vutl.nlme4)$fixed[,c(1,2)], 1, diff)
fixASi <- c(fixef(Vutl.nlme4)[1], fixef(Vutl.nlme4)[1] + fixef(Vutl.nlme4)[c(2:4)])
cbind(fixASi, fixASi + cis[c(1:4)] %*% cbind(-1, 1))
fixASf <- c(fixef(Vutl.nlme4)[5], fixef(Vutl.nlme4)[5] + fixef(Vutl.nlme4)[c(6:8)])
cbind(fixASf, fixASf + cis[c(5:8)] %*% cbind(-1, 1))
fixtaud <- c(fixef(Vutl.nlme4)[9], fixef(Vutl.nlme4)[9] + fixef(Vutl.nlme4)[c(10:12)])
cbind(fixtaud, fixtaud + cis[c(9:12)] %*% cbind(-1, 1))

AS.fixed <- rbind(	cbind(fixASi, fixASi + cis[c(1:4)] %*% cbind(-1, 1)),
										cbind(fixASf, fixASf + cis[c(5:8)] %*% cbind(-1, 1)),
										cbind(fixtaud, fixtaud + cis[c(9:12)] %*% cbind(-1, 1)),
											cbind(fixef(Vutl.nlme4)[13], fixef(Vutl.nlme4)[13] + cis[c(13)] %*% cbind(-1, 1))
											)
AS.fixed <- as.data.frame(AS.fixed)		
names(AS.fixed)<-c("Est","lower","upper")
AS.fixed

save(Vutl.nlme4, pred.Vutl.nlme4, AS.fixed, Vutl.noAdb, file="082005NaturePlantsnlmeActivationState.Rdata")

#################################################################################
#SETTING ASIDE THE BELOW FOR NOW
##################################################
##building 95% CI for this model based on bootstrap
###################################################
#
##code to produce bootstrapped confidence intervals for time series
##based on https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme
##sampling with replacement from among the Replicate.ST
##then sampling with replacement for the residuals in that set
##the model will subsequently be re-fit to these datasets and re-predicted
##and the set of predictions used to generate the CI!
###################################
##for each replicate (single leaf timeseries):
##the predictions and residuals are obtained
##the residuals are sampled with replacement to produce a bootstrapped set of residuals
##the boostrapped residuals are added to the predictions
##the sets of bootstrapped values for individual replicates
##(a resampled dataset with the same error characteristics as the original)
##are bound to produce a single set
##this represents one resampling of the model fit
#
##Because the original timesequences are not 100% overlapping between replicates,
##rather than use predict, the method here uses fitted/residuals
##aligned against a common timeseries backbone pulled from a square dataset
#
#hds<-c("TimePoint.s", "s.0", "s.d", "s.a", "Replicate.ST", "Genotype")
#
#######################################################
##think I need to revisit this to handle genotype level...
#resamp <- function(mod, level, dataset, N){
#
#	ds <- dataset
#	ds$Replicate.ST <- factor(as.character(ds$Replicate.ST))
#	row.names(ds) <- c(1:nrow(ds))
#
#	ts.unique <- unique(ds$TimePoint.s)
#	ts.unique <- ts.unique[order(ts.unique)]
#	
#	ts.sqr <- ds$TimePoint.s
#	ts.samp <- ds$TimePoint * 60#this is converting the original unit values, which include NAs, to s
#	
#	fit <- fitted(mod, level = level)
#	fit.sqr <- replace(ts.samp, !is.na(ts.samp), fit)
#	
#	res <- resid(mod, level = level)
#	res.sqr <- replace(ts.samp, !is.na(ts.samp), res)
#	
#	ts.sqr.byrep <- by(ts.sqr, ds$Replicate.ST, identity)
#	fit.sqr.byrep <- by(fit.sqr, ds$Replicate.ST, identity)
#	res.sqr.byrep <- by(res.sqr, ds$Replicate.ST, identity)
#	
#	#establish a sequence of random draws with replacement that sample N times from within the set of replicates 
#	r.set <- replicate(N, sample.int(length(res.sqr.byrep), length(res.sqr.byrep), replace = TRUE))
#
#	#function that pulls resids, resamples them with replacement
#	#adds them back to their respective fitted vals
#	#then appends them to a dataframe
#	#. is a resampled index for the list by replicate
#	shuffle <- function(.){
#		rs.rep <- res.sqr.byrep[.]
#		#function to resample non-NA values into 
#		rs.rep.nona <- sapply( rs.rep, function(.){ .[!is.na(.)] } )
#		rs.res.nona <- unlist(sapply( rs.rep.nona, function(.){ sample(., size = length(.), replace = TRUE) } ))
#		rs.sqr <- unlist(rs.rep)
#		rs.sqr <- replace(rs.sqr, !is.na(rs.sqr), rs.res.nona)
#		rs.df <- data.frame(ds[, hds], ActivationState = fit.sqr+rs.sqr)
#		rs.df[order(rs.df$TimePoint.s, rs.df$Replicate.ST), ]
#	}
#	
#	apply(r.set, 2, shuffle)
#		
#}
#
##test
#resamp(Vutl.nlme4, 1, Vutl.noAdb, 1)
#
#plme <- function(., mod, level){
#
#	u.me <- try(
#							update(mod,
#											data = .,
#											start = list( fixed = fixef(mod), random = as.matrix(ranef(mod)) )
#											)
#							)
#
#	if (attr(u.me, "class") != "try-error") {	pp <- predict(u.me, level = 0)
#											replace(.$ActivationState, !is.na(.$ActivationState), pp)
#											} else { replace(.$ActivationState, !is.na(.$ActivationState), NA) }
#
#}
#
##function to pull out a quantiles from output of resamp
#get.CI <- function(., probs = c(0.025, 0.975)){
#	ci <- data.frame(t(apply(., 1, quantile, probs, na.rm = TRUE)))
#	names(ci) <- c("lower", "upper")
#	ci
#}
#
#sq.it <- function(x, y){
#	replace(x$TimePoint, !is.na(x$TimePoint), y)
#	}
#
#p.it <- function(x, y, level){
#	predict(x, newdata = y[, hds], level = level)
#	}
#	
#p.it.df <- function(x, y, level){
#	data.frame(pred = predict(x,
#														newdata = y[, hds],
#														level = level
#														),
#							y[, hds]
#							)
#	}
#
##use the above to append 95%CI for the mean prediction to the dataset
##the CI is generated using 1000 replicates
##the 500 reps used in the previous code worked for IT82E16, but that model fails less than some of the others when bootstrapping...
##so 100 used aiming to get N~400 per genotype
#
##group level predicted values at TimePoints.s
#p.rep <- p.it(Vutl.nlme4, Vutl.noAdb, level = 1)
##order to match output from resamp
#p.rep <- p.rep[order(Vutl.noAdb$TimePoint.s, Vutl.noAdb$Replicate.ST)]
##names(p.rep) <- paste0("byrep_", names(p.rep))
#
##a fully smooth set of predictions for plotting
#butry.TP.s <- c((-4 * 60):(40 * 60))
#butry.s.0 <- ifelse(butry.TP.s < 0, 1, 0)
#butry.s.d <- ifelse(butry.TP.s < 0 | butry.TP.s >=  1200, 0, 1)
#butry.s.a <- ifelse(butry.TP.s >=  1200, 1, 0)
#butry.R.ST <- unique(Vutl$Replicate.ST)
#butry.newdata <- data.frame(
#														TimePoint.s = rep(butry.TP.s, length(butry.R.ST)), 
#														s.0 = rep(butry.s.0, length(butry.R.ST)), 
#														s.d = rep(butry.s.d, length(butry.R.ST)), 
#														s.a = rep(butry.s.a, length(butry.R.ST)), 
#														Replicate.ST = rep(butry.R.ST, each = length(butry.TP.s))
#														)
#butry.newdata$Genotype <- sapply(strsplit(as.character(butry.newdata$Replicate.ST), split = "_"), function(.){ .[1] } )
#butry.newdata$Genotype <- factor(butry.newdata$Genotype, levels = levels(Vutl$Genotype))
##note levels needs to be specified on above line or the next function will re-factor this with names in wrong order
##list.butry.newdata <- by(butry.newdata, butry.geno, identity)
#
#butry.p.rep <- p.it.df(Vutl.nlme4, Vutl.noAdb, level = 1)
##names(butry.p.rep) <- paste0("smooth_", names(p.rep))
#
##fixef level predicted values at TimePoints.s
#p.geno <- p.it(Vutl.nlme4, Vutl.noAdb, level = 0)
##order to match output from resamp
#p.geno <- p.geno[order(Vutl.noAdb$TimePoint.s, Vutl.noAdb$Replicate.ST)]
##names(p.geno) <- paste0("fitted_", names(p.geno))
#
#butry.p.geno <- p.it.df(Vutl.nlme4, butry.newdata, level = 0)
##names(butry.p.geno) <- paste0("smooth_", names(p.geno))
#
##resample
#set.seed(100)
#rs <- resamp(Vutl.nlme4, Vutl.noAdb, level = 1, N = 1000)
##fits nlme models... very extended...
#p.rs <- lapply(rs, plme, mod = Vutl.nlme4, level = 0)#note that these models fail in more than 50% of cases
#p.rs.format <- do.call(cbind, p.rs)
##check number of bootstrap reps achieved - should be > = ~400
##very rough estimate because NAs are present in some models anyway
#apply(p.rs.format, 2, function(.){ all(is.na(.)) } )#998?
#CI95 <- t(apply(p.rs.format, 1, quantile,  probs = c(0.025, 0.975), na.rm=TRUE))
##make a dataframe with raw data and predictions (both of which are ordered according to the original dataframes, not by TimePoint.s as the CI95 are
##p.rep.sq <- mapply(sq.it, list.data, p.rep)
##p.geno.sq <- mapply(sq.it, list.data, p.geno)
#
#data.frame(p.geno,p.rep,CI95)#################somethings's not right...possibly missing geno term in predict for
#
#
#list.preds <- mapply(data.frame, list.data, p.geno = p.geno, p.rep = p.rep, SIMPLIFY  =  FALSE, USE.NAMES  =  TRUE)
##reorder this to match CI95
#list.preds <- lapply(list.preds, function(.){.[order(.$TimePoint.s, .$Replicate.ST), ]})
##now merge with CI95
#list.preds <- mapply(data.frame, list.preds, list.CI95, SIMPLIFY  =  FALSE, USE.NAMES  =  TRUE)
##and get rid of NA for plotting purposes
#list.preds.omit <- lapply(list.preds, na.omit)
#
#
##save.image("072013timesequencemodel_gasexchange.Rdata")