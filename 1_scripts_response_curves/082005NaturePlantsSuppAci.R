#here, re-drawing the A/ci responses, just for those curves that were needed to do one-point Vcmax

rm(list = ls())

setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rdata")
setwd("~/Downloads")

load("082005NaturePlantsnlmeVcmax.Rdata")
load("082005summarize.Rdata")

str(AC.fits_GammastarKcogm_fixed_paired)
names(AC.fits_GammastarKcogm_fixed_paired)

head(ilGKg.df)
unique(ilGKg.df$plant)#matches levels

plot.ACfit2 <- function(input.AC){

	plot.Acifit(input.AC, xlim = c(0, 90), ylim = c(-5, 55), PhiPSII = TRUE, F_Slimitation = TRUE)
	title(main = input.AC$data$plant[1])
	#plot.Acifit(input.AC, xlim = c(0, 20), ylim = c(-5, 20), PhiPSII = FALSE, F_Slimitation = FALSE)
	#title(main = input.AC$data$plant[1])

	}

plot.ACfit.bygeno <- function(AC.fits){
	#line modified compared to 082005summarize.Rdata
	for (i in levels(ilGKg.df$geno)){
	#removed nar specification here as it doesn't matter with a single fit
		par(mfrow = c(3, 2), las = 1)
		lapply(AC.fits[grep(i, names(AC.fits), value = TRUE)], plot.ACfit2)
	}
}

inboth <- names(AC.fits_GammastarKcogm_fixed_paired) %in% unique(ilGKg.df$plant)
inboth <- names(AC.fits_GammastarKcogm_fixed_paired)[inboth]
inboth <- AC.fits_GammastarKcogm_fixed_paired[inboth]

pdf("082005NaturePlantsSuppAci.pdf", w = 6.5, h = 8)
plot.ACfit.bygeno(inboth)
dev.off()
