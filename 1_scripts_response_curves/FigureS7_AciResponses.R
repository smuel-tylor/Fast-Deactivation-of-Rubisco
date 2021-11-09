################################################################################
#Re-drawing the A/ci responses
# just for those curves that were used to do one-point Vcmax

library(here)

load(here("output/nlmeVcmax.Rdata"))
load(here("output/summarize.Rdata")) #includes plot.Acifit function

#plotting functions
plot.ACfit2 <- function(input.AC){
  
  plot.Acifit(input.AC,
              xlim = c(0, 90),
              ylim = c(-5, 53),
              PhiPSII = TRUE,
              F_Slimitation = FALSE
  )
  points(input.AC$F_Slimitation["Pci.op"],
         input.AC$F_Slimitation["A.op"],
         col = "red",
         pch = 19
         )
  title(main = input.AC$data$plant[1])
}

plot.ACfit.bygeno <- function(AC.fits){
	#line modified compared to summarize.Rdata
  for (i in unique(
    sapply(strsplit(names(ACfits), "_"), function(.){.[1]} )
  )
  ){
	#removed mar specification here as it doesn't matter with a single fit
		par(mfrow = c(3, 2),
		    las = 1,
		    tcl = 0.3,
		    mgp = c(2, 0.5, 0)
		    )
		lapply(AC.fits[grep(i, names(AC.fits), value = TRUE)], plot.ACfit2)
	}
}

#subset and reformat data
names(AC.fits_GammastarKcogm_fixed_paired)
#levels for plants used in nlmeVcmax.R
unique(ilGKg.df.sub$plant)

inboth <- names(AC.fits_GammastarKcogm_fixed_paired) %in%
  unique(ilGKg.df.sub$plant)

ACfits <- AC.fits_GammastarKcogm_fixed_paired[inboth]
names(ACfits) <- gsub("TVNu-1948", "V. sp. Savi", names(ACfits))

#function to replace plant
rep_plcode <- function(input.AC){
  input.AC$data$plant <- gsub("TVNu-1948", "V. sp. Savi", input.AC$data$plant)
  #for tidiness
  input.AC$data$geno <- gsub("TVNu-1948", "V. sp. Savi", input.AC$data$geno)
  input.AC
  }
ACfits <- lapply(ACfits, rep_plcode)  

pdf(here("output/FigureS7.pdf"),
    w = 16 / 2.58, h = 21.5 / 2.58
    )
plot.ACfit.bygeno(ACfits)
dev.off()
