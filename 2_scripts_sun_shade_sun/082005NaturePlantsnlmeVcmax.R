#Non-linear-mixed-effects model building for the coefficients needed in diurnal model
#half times for Vcmax and ActivationState, and light response curves

#The order of the "Genotype" factor is adjusted to be: Adenantha,  TVNu-1948,  IT82E-16,  IT86D-1010
#and square datasets produced for subsequent use in bootstrapping CIs

#clear workspace
rm(list = ls())

library(lattice)
library(nlme)

setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rdata")
setwd("~/Downloads")
load("082005summarize.Rdata")
objects()

#to make the piecewise model work,  I need to add appropriate time sequence factors to the dataset

#work through the modelling with "inds.list_GammastarKcogm_fixed_paired"
#the best fitting A/ci analysis
 
#produce a mirror object with a simpler names
ilGKg <- inds.list_GammastarKcoInfgm_fixed_paired
str(ilGKg)

#The only effective way to set limits on the Vcmax controlled portion of the biochemistry
#is A < Aj | Ci

AC.j <- function(ACisumm, ind){
	
	J = ACisumm$chosen.mod$J
	Gamma.star = ACisumm$chosen.mod$Gamma.star
	Rd = ACisumm$chosen.mod$Rd
	gm = ACisumm$chosen.mod$gm
	Pci.J = ind$Pci
	stoich = "NADPH"
					
	ind$Aj.Pci <- RuBP.limited(J, Gamma.star, Rd, gm, Pci.J, stoich)
	#ind$Vcmax.f <- ACsumm$chosen.mod$Vcmax#adding this because I decided I need it later to set asymptote of induction model
	
	ind
}

ilGKg <- mapply(AC.j, AC.fits_GammastarKcogm_fixed_paired, ilGKg, SIMPLIFY = FALSE)

#simplify this to a dataframe, 
ilGKg.df <- data.frame(do.call(rbind, ilGKg))
#recode the timings to be zero at the start of the shade period
ilGKg.df$induction.s <- ilGKg.df$induction.s + 1200

pdf("082005NaturePlantsnlmeVcmax_assessJlimitation.pdf", paper = "a4", w = 8, h = 11)
header.page(main = expression(
															atop("Does "*italic(A)*"(black line) intersect with",
																		italic(A)[J]*" (red line),  predicted at steady state "*italic(J)*" & "*italic(c)[i]*" during induction")
																		)
															)

for (i in unique(ilGKg.df$geno)){
	par(mfrow = c(3, 2), mar = c(4, 5.5, 2, 0.5), las = 1, cex = 1.2)
		lapply(ilGKg[grep(i, names(ilGKg),
						value = TRUE)],
						function(.){ plot(A~induction.s, data=., 
															xlab = "time since end of shade (s)", 
															ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1)), 
															xlim = c(-120, 900), ylim = c(0, 40), 
															type = "l", main = .$plant[1])
												lines(Aj.Pci ~ induction.s, data = ., col = 2)
												lines(rep(60, 2), c(-1e6, 1e6), col = 5)
												lines(rep(300, 2), c(-1e6, 1e6), col = 3)
												lines(rep(600, 2), c(-1e6, 1e6), col = 4)
												}
			)
}

dev.off()

#don't know that 5 min will be enough for good fits,  but if it works,  it's the most easily justified...
#probably also need to look at 10 min...

#limit the timesequences to keep them consistent
ilGKg.df <- ilGKg.df[ilGKg.df$induction.s >= -120&ilGKg.df$induction.s <= 1200 + 300, ]
#and get rid of genotypes not used for the Rubisco activation state work
levels(ilGKg.df$geno)
levels(ilGKg.df$geno)[c(5, 6)] <- NA
ilGKg.df <- ilGKg.df[!is.na(ilGKg.df$geno), ]

#retain only the columns actually needed for the model - necessary to ensure an na.omit below works OK because FLR values are usually all NA
ilGKg.df <- ilGKg.df[, c("induction.s", "geno", "plant", "Astar", "Vcmax.t")]#, "Vcmax.f")]

#Make the dataset square so that resampling for CI will work
#pad out the missing values in each set of TimePoints*Replicate.ST
#Vutl$Replicate.ST.s <- paste(Vutl$Replicate.ST, Vutl$TimePoint.s, sep = "_")
#a vector that gives what would be expected in a square dataset
Ts <- unique(ilGKg.df$induction.s)
pl <- unique(ilGKg.df$plant)
E.plTs <- paste(rep(pl, each = length(Ts)), rep(Ts, times = length(pl)), sep = "_")

ilGKg.df2 <- data.frame(geno = sapply(strsplit(E.plTs, "_"), function(.){ .[1] }), 
												plant = sapply(strsplit(E.plTs, "_"), function(.){ paste(.[c(1, 2)], collapse = "_") }), 
												induction.s = sapply(strsplit(E.plTs, "_"), function(.){ .[3] })
												)
head(ilGKg.df2)

nrow(ilGKg.df)
nrow(ilGKg.df2)
nrow(merge(ilGKg.df2, ilGKg.df, all.x = TRUE))#seems to work nicely
table(merge(ilGKg.df2, ilGKg.df, all.x = TRUE)$plant)#as above

ilGKg.df <- merge(ilGKg.df2, ilGKg.df, all.x = TRUE)
#clean up
rm(ilGKg.df2, Ts, pl, E.plTs)

#fix all the relevant factors by exporting this dataset and re-reading
#not actually different from previously dated versions of the same file, but changing the name here to avoid confusion
write.csv(ilGKg.df, "Vu_timelapse_gasexchange_082005_square.csv",row.names=FALSE)
ilGKg.df <- read.csv("Vu_timelapse_gasexchange_082005_square.csv")

#update genotype order
ilGKg.df$geno <- factor(ilGKg.df$geno, levels = c("Vadenantha",  "TVNu-1948",  "IT82E-16",  "IT86D-1010"))
#to make everything run smoothly,  the dataframe now needs to be reordered to match this!
ilGKg.df <- ilGKg.df[order(ilGKg.df$geno), ]

#add a 'time' variable that includes NA values where samples were not collected (used below)
ilGKg.df$inductionNA.s <- replace(ilGKg.df$induction.s, is.na(ilGKg.df$Vcmax.t), NA)

xyplot(Astar*5 + Vcmax.t ~ induction.s | plant, data = ilGKg.df)
xyplot(Astar ~ Vcmax.t | plant, data = ilGKg.df)
#Makes sense,  since Astar is assumed to be Vcmax limited and all other parameters are held constant
#the two are directly proportional

#Modelling Vcmax.t is going to be more obviously parallel with Rubisco activation state
#Astar is equivalent to Vcmax.t,  and is only correct over the Vcmax limited range
#(true Astar will be smaller in J limited range)  

#####################
#timeseries analyses#
#####################
#In a previous iteration, I realised this is only valid for the initial post-shade phase
#add switching coefficients
#noting that unlike the Rubisco activation state experiment,  I cannot model the 0-1200 period
#ilGKg.df$s.0 <- ifelse(ilGKg.df$induction.s < 0,  1,  0)
#Vutl$s.d <- ifelse(Vutl$TimePoint.s < 0 | Vutl$TimePoint.s >= 1200,  0,  1)
ilGKg.df$s.a <- ifelse(ilGKg.df$induction.s > 1260,  1,  0)#fails if >=1260!
ilGKg.df.ind <- ilGKg.df[ilGKg.df$s.a == 1, ]

ind.nlsList <- nlsList(Vcmax.t ~ (Vcmax.f - (Vcmax.f - Vcmax.i) * exp( -(induction.s - 1200) / tau.a)) | plant, 
												start = c(Vcmax.i = 60, Vcmax.f = 170, tau.a = 120), 
												data = ilGKg.df.ind, na.action = na.exclude)
summary(ind.nlsList)
plot(ind.nlsList,  plant ~ resid(.),  abline  =  0 )
plot(intervals(ind.nlsList))
#some really poorly constrained estimates of the coefs 
ind.nlsList2 <- update(ind.nlsList,
												data = ilGKg.df[ilGKg.df$plant!="IT86D-1010_6"&ilGKg.df$plant!="IT82E-16_2", ]
												)
plot(intervals(ind.nlsList2))

#presumably though,  some level of pooling is needed to keep residuals homoskedasitc and normal
ind.nls <- nls(Vcmax.t ~ (Vcmax.f - (Vcmax.f - Vcmax.i) * exp(-(induction.s - 1200) / tau.a)), 
								start = c(Vcmax.i = 60, Vcmax.f = 170, tau.a = 120), 
								data = ilGKg.df.ind, na.action = na.exclude
								)
summary(ind.nls)
#residual is much larger than for the fully parameterised model
plot(ind.nls,  plant ~ resid(.),  abline  =  0 )
#explains why

ind.nlme <- nlme(ind.nlsList)
ind.nlme
#too many random effects...
plot(ranef(ind.nlme))
#but no obvious way to structure them
#so,  the next thing to attempt would be fixed effects by genotype
ind.nlme2 <- update(ind.nlme,
											fixed = Vcmax.i + Vcmax.f + tau.a ~ geno,
											start = c(fixef(ind.nlme)[1], 0, 0, 0,
																fixef(ind.nlme)[2], 0, 0, 0,
																fixef(ind.nlme)[3], 0, 0, 0
																),
											control = lmeControl(maxIter=100)
										)
anova(ind.nlme, ind.nlme2)
#this is a bad step... AIC increases
anova(ind.nlme2)#possibly only Vcmax.i is important as a fixed effect?
#or,  drop tau.a as a genotype-level effect because it clearly isn't differentiating
ind.nlme3 <- update(ind.nlme,
										fixed = list(Vcmax.i + Vcmax.f ~ geno, tau.a~1),
										start = c(fixef(ind.nlme)[1], 0, 0, 0,
															fixef(ind.nlme)[2], 0, 0, 0,
															fixef(ind.nlme)[3]
															),
										control = lmeControl(maxIter=100)
										)
anova(ind.nlme, ind.nlme3)#smaller but also bad shift in AIC,  similar effect on likelihood
pairs(ind.nlme3)
anova(ind.nlme3)#looks like all fixef stay
ind.nlme4 <- update(ind.nlme,
										fixed = list(Vcmax.i + Vcmax.f ~ geno, tau.a~1),
										start = c(fixef(ind.nlme)[1], 0, 0, 0,
															fixef(ind.nlme)[2], 0, 0, 0,
															fixef(ind.nlme)[3]
															),
										random = Vcmax.i + Vcmax.f ~ 1
										)
ind.nlme5 <- update(ind.nlme,
										fixed = list(Vcmax.i + Vcmax.f ~ geno, tau.a~1),
										start = c(fixef(ind.nlme)[1], 0, 0, 0,
															fixef(ind.nlme)[2], 0, 0, 0,
															fixef(ind.nlme)[3]
															), 
										random = Vcmax.i + tau.a ~ 1
										)
anova(ind.nlme3, ind.nlme5)
anova(ind.nlme3, ind.nlme4)
#better off dropping tau.a as a random effect					
anova(ind.nlme4)#now geno is not differentiating Vcmax
ind.nlme6 <- update(ind.nlme,
										fixed = list(Vcmax.i ~ geno, Vcmax.f + tau.a ~ 1),
										start = c(fixef(ind.nlme)[1], 0, 0, 0,
															fixef(ind.nlme)[2],
															fixef(ind.nlme)[3]
															), 
										random = Vcmax.i + Vcmax.f ~ 1
										)
anova(ind.nlme4, ind.nlme6)
anova(ind.nlme6)
intervals(ind.nlme6)
plot(augPred(ind.nlme6, primary = ~induction.s, level = 0:1))#this is nice

#produce a nicer plot of this to match ActivationState
t.smth <- c(1200:1800)
p.smth <- as.character(unique(ilGKg.df.ind[ , "plant"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.Vc <- data.frame(induction.s = rep(t.smth, length(p.smth)), 
											geno = rep(g.smth, each = length(t.smth)), 
											plant = rep(p.smth, each = length(t.smth))
											)
nd.Vc$plant <- factor(nd.Vc$plant, levels = levels(ilGKg.df.ind$plant))
nd.Vc$geno <- factor(nd.Vc$geno, levels = levels(ilGKg.df.ind$geno))

pred.ind.nlme6 <- data.frame(Vc.geno = predict(ind.nlme6, newdata = nd.Vc, level = 0), 
															Vc.plant = predict(ind.nlme6, newdata = nd.Vc, level = 1), 
															plant = nd.Vc$plant,
															geno = nd.Vc$geno, 
															induction.s = nd.Vc$induction.s 
															)
															
fhds <- c("plant", "induction.s", "Vcmax.t")
facs.ind.nlme6 <- ilGKg.df.ind[, fhds]
pred.ind.nlme6 <- merge(pred.ind.nlme6, facs.ind.nlme6, all.x = TRUE)

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
				ylim = c(0, 260),
				xlab = "time from end of shade (s)",
				ylab = expression(italic(V)["c,max"]~~(mu*mol~~m^-2~~s^-1)),
				main = ordp$plant[1],
				type = "n",
				axes = FALSE)
	axis(side = 1, at = seq(1200, 1800, 120), labels = seq(1200, 1800, 120) - 1200, las = 1)
	axis(side = 2, at = seq(0, 250, 50), las = 1)
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

pdf("082005NaturePlantsnlmeVcmax_plot_nlme6.pdf", w = 6, h = 9)							

par(mar = c(4.5, 5.5, 3, 1), tcl = 0.4, oma = c(0, 0, 0, 0), cex = 1.2, cex.lab = 1.2, cex.axis = 1.2)

glist <- levels(pred.ind.nlme6$geno)
lapply(glist, allrep_plot, p.all = pred.ind.nlme6)

dev.off()

#Here's a plot of the initial part of the Vcmax ~ t response
#showing that there are only very marginal effects after the initial 5 min period when Vcmax is always limiting
pdf("082505NaturePlantsnlmeVcmax_extrapolation.pdf", paper = "a4", w = 8, h = 8)
nd <- data.frame(induction.s = rep(seq(1200, 2400, 10), 4),
									geno = rep(levels(ilGKg.df$geno), each = length(seq(1200, 2400, 10)))
									)
in6g <- predict(ind.nlme6, 
								newdata = nd, 
								level = 0
								)
in6g <- data.frame(nd, Vcmax.t = in6g)
par(mfrow=c(1,1),mar = c(5, 5, 5, 5))
plot(Vcmax.t ~ induction.s,
			data = ilGKg.df[ilGKg.df$induction.s > 1260, ], 
			las = 1,
			pch = 21,
			col = NULL,
			bg = rgb(0, 0, 0, alpha = 0.2), 
			xlab = "time since shade (min)",
			xlim = c(1200, 1800), xaxt = "n", 
			ylab = expression(italic(V)["c, max"]~~(mu*mol~~m^-2~~s^-1)),
			ylim = c(0, 225)
			)
axis(side = 1,
			at = c(1200, 1500, 1800),
			labels = c(0, 5, 10)
			)
title(main = expression(
												atop(italic(V)["c, max"]*" timeseries from gas exchange in four "*italic(Vigna)*" genotypes",
															"effect of projection beyond "*italic(V)["c.max"]*" limited range"
															)
											  )
			)

by(in6g,
		in6g$geno,
		function(.){ lines(Vcmax.t ~ induction.s, data = .)
									points(30 * diff(.$Vcmax.t) ~ .$induction.s[c(2:nrow(.))],
													col = rgb(1, 0, 0, alpha = 0.3),
													pch = 21
													)
									}
		)
#points(30*diff(177-(177-88)*exp(-seq(0, 600, 10)/129))~seq(1210, 1800, 10), col = 2, pch = 19)
axis(side = 4,
			at = c(0, 60, 120, 180),
			labels = c(0, 60, 120, 180) / 30,
			col = 2,
			col.axis = 2,
			las = 1
			)
mtext(side = 4,
			line = 3,
			expression("change in "*italic(V)["c, max"]*" per 10 s interval"),
			col = 2
			)
			
dev.off()


#fixed effects CIs as one-tailed value
cis <- apply(intervals(ind.nlme6)$fixed[,c(1,2)], 1, diff)
fixVcmax <- c(fixef(ind.nlme6)[1], fixef(ind.nlme6)[1] + fixef(ind.nlme6)[c(2:4)])
cbind(fixVcmax, fixVcmax + cis[c(1:4)] %*% cbind(-1, 1))

Vcind.fixed <- rbind(cbind(fixVcmax, fixVcmax + cis[c(1:4)] %*% cbind(-1, 1)),
											cbind(fixef(ind.nlme6)[5:6], fixef(ind.nlme6)[5:6] + cis[c(5:6)] %*% cbind(-1, 1))
											)
Vcind.fixed <- as.data.frame(Vcind.fixed)		
names(Vcind.fixed)<-c("Est","lower","upper")
Vcind.fixed


save(in6g,ind.nlme6,Vcind.fixed,ilGKg.df,file="082005NaturePlantsnlmeVcmax.Rdata")

##########################################
##Setting aside the below for now
########################################################
##building 95% CI for this model based on bootstrapping
########################################################
#
##code to produce bootstrapped confidence intervals for time series
##based on https://stats.stackexchange.com/questions/231074/confidence-intervals-on-predictions-for-a-non-linear-mixed-model-nlme
##sampling with replacement from among the Replicate.ST
##then sampling with replacement for the residuals in that set
##the model will subsequently be re-fit to these datasets and re-predicted
##and the set of predictions used to generate the CI!
###################################
##for each of replicate:
##the predictions and residuals are obtained
##the residuals are sampled with replacement to produce a bootstrapped set of residuals
##the boostrapped residuals are added to the predictions
##the four sets of bootstrapped values
##(representing a resampled dataset with the same error characteristics as the original)
##are bound to produce a single set
##this therefore represents one resampling of the model fit
#
##Because the original timesequences are not 100% overlapping between replicates,
##rather than use predict, the method here uses fitted/residuals
##aligned against a common timeseries backbone pulled from a square dataset
#resamp<-function(mod,level,dataset,N){
#
#	ds<-dataset
#	ds$Replicate.ST<-factor(as.character(ds$Replicate.ST))
#	row.names(ds)<-c(1:nrow(ds))
#
#	ts.unique<-unique(ds$TimePoint.s)
#	ts.unique<-ts.unique[order(ts.unique)]
#	
#	ts.sqr<-ds$TimePoint.s
#	ts.samp<-ds$TimePoint*60#this is converting the original unit values, which include NAs, to s
#	
#	fit<-fitted(mod,level=level)
#	fit.sqr<-replace(ts.samp,!is.na(ts.samp),fit)
#	
#	res<-resid(mod,level=level)
#	res.sqr<-replace(ts.samp,!is.na(ts.samp),res)
#	
#	ts.sqr.byrep<-by(ts.sqr,ds$Replicate.ST,identity)
#	fit.sqr.byrep<-by(fit.sqr,ds$Replicate.ST,identity)
#	res.sqr.byrep<-by(res.sqr,ds$Replicate.ST,identity)
#	
#	#establish a sequence of random draws with replacement that sample N times from within the set of replicates 
#	r.set<-replicate(N,sample.int(length(res.sqr.byrep),length(res.sqr.byrep),replace=TRUE))
#
#	#function that pulls resids, resamples them with replacement
#	#adds them back to their respective fitted vals
#	#then appends them to a dataframe that can be used for
#	#. is a resampled index for the list by replicate
#	shuffle<-function(.){
#		rs.rep<-res.sqr.byrep[.]
#		#function to resample non-NA values into 
#		rs.rep.nona<-sapply( rs.rep, function(.){ .[!is.na(.)] } )
#		rs.res.nona<-unlist(sapply( rs.rep.nona, function(.){ sample(., size=length(.), replace=TRUE) } ))
#		rs.sqr<-unlist(rs.rep)
#		rs.sqr<-replace(rs.sqr,!is.na(rs.sqr),rs.res.nona)
#		rs.df<-data.frame(ds[,c("TimePoint.s","s.0","s.d","s.a","Replicate.ST")],ActivationState=fit.sqr+rs.sqr)
#		rs.df[order(rs.df$TimePoint.s,rs.df$Replicate.ST),]
#	}
#	
#	apply(r.set,2,shuffle)
#		
#}
#
##test
#resamp(nlts[[1]],1,list.data[[1]],1)
#
#plme<-function(.,mod,level){
#
#	u.me<-try(update(mod,data=.,start=list(fixed=fixef(mod),random=as.matrix(ranef(mod)))))
#
#	if (attr(u.me,"class")!="try-error") {	pp<-predict(u.me,level=0)
#											replace(.$ActivationState,!is.na(.$ActivationState),pp)
#											} else { replace(.$ActivationState,!is.na(.$ActivationState),NA) }
#
#}
#
##function to pull out a quantiles from output of resamp
#get.CI<-function(.,probs=c(0.025,0.975)){
#	ci<-data.frame(t(apply(.,1,quantile,probs,na.rm=TRUE)))
#	names(ci)<-c("lower","upper")
#	ci
#}
#
#sq.it<-function(x,y){
#	replace(x$TimePoint, !is.na(x$TimePoint), y)
#	}
#
#p.it<-function(x,y,level){
#	predict(x,newdata=y[,c("TimePoint.s","s.0","s.d","s.a","Replicate.ST")],level=level)
#	}
#	
#p.it.df<-function(x,y,level){
#	data.frame(pred=predict(x,newdata=y[,c("TimePoint.s","s.0","s.d","s.a","Replicate.ST")],level=level),y[,c("TimePoint.s","s.0","s.d","s.a","Replicate.ST")])
#	}
#
##use the above to append 95%CI for the mean prediction to the dataset
##the CI is generated using 1000 replicates
##the 500 reps used in the previous code worked for IT82E16, but that model fails less than some of the others when bootstrapping...
##so 100 used aiming to get N~400 per genotype
#
##group level predicted values at TimePoints.s
#p.rep<-mapply(p.it,nlts,list.data,MoreArgs=list(level=1))
#names(p.rep)<-paste0("byrep_",names(p.rep))
#
##a fully smooth version for plotting
#butry.TP.s<-c((-4*60):(40*60))
#butry.s.0<-ifelse(butry.TP.s < 0, 1, 0)
#butry.s.d<-ifelse(butry.TP.s < 0 | butry.TP.s >= 1200, 0, 1)
#butry.s.a<-ifelse(butry.TP.s >= 1200, 1, 0)
#butry.R.ST<-unique(Vutl$Replicate.ST)
#butry.newdata<-data.frame(
#							TimePoint.s=rep(butry.TP.s,length(butry.R.ST)),
#							s.0=rep(butry.s.0,length(butry.R.ST)),
#							s.d=rep(butry.s.d,length(butry.R.ST)),
#							s.a=rep(butry.s.a,length(butry.R.ST)),
#							Replicate.ST=rep(butry.R.ST,each=length(butry.TP.s))
#							)
#butry.Rep.ST=as.character(butry.newdata$Replicate.ST)
#butry.geno<-sapply(strsplit(butry.Rep.ST,split="_"), function(.){ .[1] } )
#butry.geno<-factor(butry.geno,levels=levels(Vutl$Genotype))#note levels needs to be specified here or the next function will re-factorise this with names in wrong order
#list.butry.newdata<-by(butry.newdata,butry.geno,identity)
#
#butry.p.rep<-mapply(p.it.df,nlts,list.butry.newdata,MoreArgs=list(level=1),SIMPLIFY=FALSE)
#names(butry.p.rep)<-paste0("smooth_",names(p.rep))
#
##fixef level predicted values at TimePoints.s
#p.geno<-mapply(p.it,nlts,list.data,MoreArgs=list(level=0))
#names(p.geno)<-paste0("fitted_",names(p.geno))
#
#butry.p.geno<-mapply(p.it.df,nlts,list.butry.newdata,MoreArgs=list(level=0),SIMPLIFY=FALSE)
#names(butry.p.geno)<-paste0("smooth_",names(p.geno))
#
##resample
#set.seed(100)
#rs<-mapply(resamp,nlts,list.data,MoreArgs=list(level=1,N=1000))
##fits nlme models... very extended...
#p.rs<-mapply(plme,rs,rep(nlts,each=1000),MoreArgs=list(level=0))#note that these models fail in more than 50% of cases
#dim(p.rs)=dim(rs)
#dimnames(p.rs)=dimnames(rs)
#p.rs.format<-lapply(dimnames(p.rs)[[2]], function(.) { do.call(cbind,p.rs[,.]) } )
##check number of bootstrap reps achieved - should be >=~400
##very rough estimate because NAs are present in some models anyway
#lapply(p.rs.format,function(.){ table(!is.na(.[1,])) } )
#list.CI95<-lapply(p.rs.format,get.CI)
##make a dataframe with raw data and predictions (both of which are ordered according to the original dataframes, not by TimePoint.s as the CI95 are
##p.rep.sq<-mapply(sq.it,list.data, p.rep)
##p.geno.sq<-mapply(sq.it,list.data, p.geno)
#
#list.preds<-mapply(data.frame,list.data,p.geno=p.geno,p.rep=p.rep,SIMPLIFY = FALSE, USE.NAMES = TRUE)
##reorder this to match CI95
#list.preds<-lapply(list.preds,function(.){.[order(.$TimePoint.s,.$Replicate.ST),]})
##now merge with CI95
#list.preds<-mapply(data.frame,list.preds,list.CI95,SIMPLIFY = FALSE, USE.NAMES = TRUE)
##and get rid of NA for plotting purposes
#list.preds.omit<-lapply(list.preds,na.omit)
#