#analysis for planned Nature Plants article
#predictions for Rubisco limitation based on gas exchange and Rubisco activation state
#an update on the RIPE Annual Meeting script in that gas exchange estimates have now been revisited and corrected

rm(list=ls())

#colour scheme
#comments from lab meeting suggest we shouldn't use the two gold colours, which are similar
#Since the RIPE palette is a bit limited in this respect, I think it makes better sense to use
#black and gray filled for the wild rels
#green and light green open for the cowpea
#I've actually tweaked the light green here so that it's gray scale distinguishable so open and closed aren't necessary
#also, indexing these directly with species names for convenience
	
o.cols <- c(
						"V. adenantha" = rgb(0.5*255, 0.5*255, 0.5*255, maxColorValue = 255), #dgry
						"V. sp. Savi." = rgb(0, 0, 0, maxColorValue = 255),#alpha=0.6*255), #black
						"IT82E-16" = rgb(135, 240, 115, maxColorValue = 255),#alpha=0.6*255), #lgrn
						"IT86D-1010" = rgb(70, 190, 110, maxColorValue = 255)#alpha=0.6*255), #dgrn
						#brwn=rgb(186,153,69,maxColorValue=255),#alpha=0.6*255),
						#gold=rgb(223,168,45,maxColorValue=255)#,alpha=0.6*255)
						)

t.cols <- c(
						"V. adenantha" = rgb(0.5*255, 0.5*255, 0.5*255, maxColorValue=255, alpha=0.2*255),
						"V. sp. Savi."  = rgb(0, 0, 0, maxColorValue=255, alpha=0.2*255),
						"IT82E-16" = rgb(135, 240, 115, maxColorValue=255, alpha=0.2*255),
						"IT86D-1010" = rgb(70, 190, 110, maxColorValue=255, alpha=0.2*255)
						)

#load up Vcmax time-series
load("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rdata/082005NaturePlantsnlmeVcmax.Rdata")

#update geno names
levels(ilGKg.df$geno)[levels(ilGKg.df$geno)=="Vadenantha"]<-"V. adenantha"
levels(ilGKg.df$geno)[levels(ilGKg.df$geno)=="TVNu-1948"]<-"V. sp. Savi."

levels(in6g$geno)[levels(in6g$geno)=="Vadenantha"]<-"V. adenantha"
levels(in6g$geno)[levels(in6g$geno)=="TVNu-1948"]<-"V. sp. Savi."

#Activation State
load("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_Biochemistry/Cowpea Time-Lapse Assays/Rmodelfitting/082005NaturePlantsnlmeActivationState.Rdata")

#update genotype names
levels(Vutl.noAdb$Genotype)[levels(Vutl.noAdb$Genotype)=="Adenantha"]<-"V. adenantha"
levels(Vutl.noAdb$Genotype)[levels(Vutl.noAdb$Genotype)=="TVNu-1948"]<-"V. sp. Savi."

levels(pred.Vutl.nlme4$Genotype)[levels(pred.Vutl.nlme4$Genotype)=="Adenantha"]<-"V. adenantha"
levels(pred.Vutl.nlme4$Genotype)[levels(pred.Vutl.nlme4$Genotype)=="TVNu-1948"]<-"V. sp. Savi."

#Nplants column width is ~8.9 cm = 3.44'
#inward ticks, sans serif

#Activation State
setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Publications/2020_Escobar et al_Rubisco_Time_Series/EGE et al Data & Figures")
pdf("EGEetal_Fig1.pdf", w = 3.5, h = 6.5)
#x11(w=3.5,h=6.5)

par(
		mfrow = c(2, 1),
		mar = c(4, 5.5, 1, 0.5),
		oma = c(0, 0, 0, 0),
		mgp = c(2, 0.5, 0),
		tcl = 0.3,
		cex.axis = 1,
		cex.lab = 1,
		las = 1
		)

plot(ActivationState ~ TimePoint.s, data = Vutl.noAdb,
			type = "n",
			ylim = c(0, 100),
			xlim = c(-5*60, 40*60),
			ylab = expression(italic(S)*" (%)"),
			xlab = expression(time~~from~~end~~of~~shade~~(min)),
			xaxt = "n"
			)

axis(side = 1, at=c(0, 10, 20, 30, 40) * 60, labels = c(-20, -10, 0, 10, 20), cex.axis = 1, lwd.ticks = 1,las = 1)
axis(side=2,at=seq(0,100,20),lwd.ticks=1,las=1)
mtext(expression(bold(a)),side=2,at=105,line=4,adj=1)

legend(1440, 20,
				xjust = 0.5,
				yjust = 0.5,
				legend = c(
										as.expression(bquote(italic(.(levels(ilGKg.df$geno)[1])))),
										bquote(italic(.(levels(ilGKg.df$geno)[2]))),
										bquote(plain(.(levels(ilGKg.df$geno)[3]))),
										bquote(plain(.(levels(ilGKg.df$geno)[4])))
										),
				pch = 21,
				pt.bg = t.cols,
				pt.lwd = 0,
				lwd = 2,
				col = o.cols,
				ncol = 1,
				bty = "n"
				)
				
lapply(	levels(Vutl.noAdb$Genotype),
				function(.){
					points(ActivationState ~ TimePoint.s,
									data = Vutl.noAdb[Vutl.noAdb$Genotype == .,],
									pch = 21,
									col = NA,
									bg = t.cols[.],
									cex = 0.7
									)
					}
)

lapply(	levels(Vutl.noAdb$Genotype),
				function(.){
					use <- pred.Vutl.nlme4[pred.Vutl.nlme4$Genotype == ., ]
					use <- use[c(1:2520), ]
					lines(AS.geno ~ TimePoint.s,
								data = use,
								col = o.cols[.],
								lwd = 2
								)
					}
)

#Activation State
plot(Vcmax.t ~ induction.s, data = ilGKg.df,
			type = "n",
			ylim = c(0, 230),
			xlim = c(20*60, 30*60),
			ylab = expression(italic(V)["c,max"]~~(mu*mol~~m^-2~~s^-1)),
			xlab = expression(time~~from~~end~~of~~shade~~(min)),
			xaxt = "n"
			)

axis(side = 1, at=c(20, 22, 24, 26, 28, 30) * 60, labels = c(0, 2, 4, 6, 8, 10), cex.axis = 1, lwd.ticks = 1,las = 1)
axis(side=2,at=seq(0,200,50),lwd.ticks=1,las=1)
mtext(expression(bold(b)),side=2,at=1.05*230,line=4.5,adj=1)

lapply(	levels(ilGKg.df$geno),
				function(.){
					points(Vcmax.t ~ induction.s,
									data = ilGKg.df[ilGKg.df$geno == . & ilGKg.df$induction.s > 1260,],
									pch = 21,
									col = NA,
									bg = t.cols[.],
									cex = 0.7
									)
					}
)

lapply(	levels(in6g$geno),
				function(.){
					use <- in6g[in6g$geno == ., ]
					lines(Vcmax.t ~ induction.s,
								data = use,
								col = o.cols[.],
								lwd = 2
								)
					}
)

dev.off()
