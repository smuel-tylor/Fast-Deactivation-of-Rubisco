#Initial cleaning and formatting of data for Nature Plants analysis

#updated 0820 for style
#updated 0121 to be consistent with use of here and renv
# with additional style tweaks in Rstudio
# and accounting for change to stringsAsFactors between R 3.x and R 4.x
# by forcing read.csv to behave as R 4.x
#Apologies, I did not use dplyr here...
#Also 0121, got rid of some old for loops and replaced with tapply
#0121, produced a new input file
#In this, I:
#- reformatted names to remove need for renaming in this script
#- removed two genotypes
#   these were not paired with biochemical analysis and therefore redundant
#script was subsequently cleaned to remove old for loops
# and unnecessary commented out lines

#for convenience plotting 
library(lattice)
library(here)

################################################################################
#load up and basic housekeeping

cp <- read.csv(here("data/RIPE_20210118_cowpeaGE.csv"),
               stringsAsFactors = FALSE
               )

#add cultivation status
cp$cult <- replace(cp$geno,
                   cp$geno == "V. adenantha" | cp$geno == "TVNu-1948",
                   "wild"
                   )
cp$cult <- replace(cp$cult, cp$cult != "wild", "cultivar")

#make factors
g.levs <- c("V. adenantha", "TVNu-1948", "IT86D-1010", "IT82E-16")

c.levs <- c("cultivar", "wild")

cp$geno <- factor(cp$geno, levels = g.levs)
cp$cult <- factor(cp$cult, levels = c.levs)

cp$block <- factor(cp$block)

#identifier for individuals
cp$plant <- factor(paste(cp$geno, cp$block, sep = "_"))
levels(cp$plant)

################################################################################
#preliminary data cleaning

################################################################################
#Aci responses
x11()
xyplot(A ~ Ci | geno, data = cp[cp$curve == "Aci", ])
x11()
xyplot(A ~ Ci, groups = geno, auto.key = T, data = cp[cp$curve == "Aci", ])
#in analysis need to identify operating points - first value in each curve
#also need to clean out non-steady state values for 430 ppm
# (might be steady state for A,  but PhiPS2 often not equivalent
# to operating point value)

#plotting to evaluate errors in dataset and identify duplicated measurements
plotAciifis <- function(dat){
  if (nrow(dat) > 1){
    dat <- dat[order(dat$Pci), ]
    plot(A ~ Pci, data = dat, main = dat$plant[1], ylim = c(-5, 55))
    points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
    axis(side = 4,
         at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4),
         labels = c(0, .1, .2, .3, .4)
    )
    mtext(expression(Phi[PSII]),
          side = 4,
          line = 3,
          at = (60 / 2) - 5,
          las = 3
    )
  } else {
    #add in a blank space where a plant was missing in a given block
    plot(1, 1,
         type = "n", axes = FALSE,
         bty = "n", main = "missing", xlab = "", ylab = "")
  }
}

#plotAciifis(AC[AC$plant == levels(AC$plant)[1], ])

plot.Acionegeno <- function(geno, data){
  x11(w = 10, h = 6)
  par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
  datg = data[data$geno == geno, ]
  lapply(levels(datg$block),
         function(.){
           plotAciifis(datg[datg$block == ., ])
           }
  )
}

AC <- cp[cp$curve == "Aci", ]

lapply(levels(AC$geno), plot.Acionegeno, data = AC)

#there are no clear issues here
#but will remove the redundant 'return to normal' points
#because these result in unequal weighting of the dataset
#around the operating point, i.e., they represent pseudoreps

#an example
AC[AC$plant == "IT82E-16_1", "co2_at"]
#in an unsorted dataset it is usually points 8 and 17

#filter out these points and plot result
AC <- by(AC, AC$plant, identity)
AC <- lapply(AC, function(.){ .[c(1:7,9:16), ] })
AC <- do.call(rbind, AC)

lapply(levels(AC$geno), plot.Acionegeno, data = AC)

################################################################################
#light responses
plotAQifis <- function(dat){
  if (nrow(dat) > 1){
    dat <- dat[order(dat$Qin), ]
    plot(A ~ Qin, data = dat, main = dat$plant[1], ylim = c(-5, 55))
    points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
    axis(side = 4,
         at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4),
         labels = c(0, .1, .2, .3, .4)
    )
    mtext(expression(Phi[PSII]),
          side = 4,
          line = 3,
          at = (60 / 2) - 5,
          las = 3
    )
  } else {
    #add in a blank space where a plant was missing in a given block
    plot(1, 1,
         type = "n", axes = FALSE,
         bty = "n", main = "missing", xlab = "", ylab = "")
  }
}

#plotAQifis(AQ[AQ$plant == levels(cp$plant)[1], ])

plot.AQonegeno <- function(geno, data){
  x11(w = 10, h = 6)
  par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
  datg = data[data$geno == geno, ]
  lapply(levels(datg$block),
         function(.){
           plotAQifis(datg[datg$block == ., ])
         }
  )
}

AQ <- cp[cp$curve == "AQ", ]

lapply(levels(AQ$geno), plot.AQonegeno, data = AQ)

#Remove 'return to 1500' points from end of curves
#an example
AQ[AQ$plant == "IT82E-16_1", "Qin"]
#in an unsorted dataset it is usually points 8 and 17

AQ <- by(AQ, AQ$plant, identity)
AQ <- lapply(AQ,
                function(.){
                  if (length(.[, 1]) > 0) {
                    .[c(1:(length(.[, 1]) - 1)), ]
                  }
                  }
                )
AQ <- do.call(rbind, AQ)

lapply(levels(AQ$geno), plot.AQonegeno, data = AQ)

################################################################################
#sun-shade-sun inductions
#look at IT86D-1010-9 and PI862537-2 inductions
# (cases where two inductions were recorded)
# - expect 2nd induction to be better
inds <- cp[cp$curve == "induction", ]

plotIndifis <- function(dat){
  if (nrow(dat) > 1){
    dat <- dat[order(dat$Qin), ]
    plot(A ~ induction.s, data = dat, main = dat$plant[1], ylim = c(-5, 55))
  } else {
    #add in a blank space where a plant was missing in a given block
    plot(1, 1,
         type = "n", axes = FALSE,
         bty = "n", main = "missing", xlab = "", ylab = "")
  }
}

#plotIndifis(inds[inds$plant == levels(inds$plant)[1], ])

plot.Indonegeno <- function(geno, data){
  x11(w = 10, h = 6)
  par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
  datg = data[data$geno == geno, ]
  lapply(levels(datg$block),
         function(.){
           plotIndifis(datg[datg$block == ., ])
         }
  )
}

lapply(levels(inds$geno), plot.Indonegeno, data = inds)

x11()
plot(gsw~induction.s,
			data = inds[inds$plant == "IT86D-1010_9", ],
			ylim = c(0, 1),
			xlim = c(-1320, 1800)
			)
plot(A~induction.s,
			data = inds[inds$plant == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
plot(Pci~induction.s,
			data = inds[inds$geno == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)

cp[cp$plant == "IT86D-1010_9" & cp$curve == "Aci", c("A", "gsw")]		

plot(A ~ Pci,
			data = AC[AC$plant == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = inds[inds$plant == "IT86D-1010_9", ],
				col = 2
				)
#Both responses overlay the A/ci response fairly well.
#The first response may better characterise induction without
# the modifications of gsw that result from Aci/AQ
#the second curve shows a more usual stomatal response during shade	
inds[inds$plant == "IT86D-1010_9" & inds$curve == "induction",
   c("obs", "induction.s")]
#drop the first induction
#i.e.,data with obs <= 363
inds <- inds[!(inds$plant == "IT86D-1010_9" & inds$obs <= 363), ]

plot(A ~ Pci,
			data = AC[AC$plant == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = inds[inds$plant == "IT86D-1010_9", ],
				col = 2
				)

plot(gsw ~ induction.s,
			data = inds[inds$plant == "IT82E-16_2", ], 
			ylim = c(0, 1),
			xlim = c(-1320, 1800)
			)
#one of these is very noisy,  but in the other gsw does not recover

plot(A ~ induction.s,
			data = inds[inds$plant == "IT82E-16_2", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
plot(Pci ~ induction.s,
			data = inds[inds$plant == "IT82E-16_2", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
#it looks like the first attempt was better despite noisy ci
#this is because gsw recovers,  whereas
#second attempt does not show complete recovery
plot(A ~ Pci,
			data = AC[cp$plant == "IT82E-16_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = inds[inds$plant == "IT82E-16_9", ],
				col = 2
				)

inds[inds$plant == "IT82E-16_2", c("obs", "induction.s", "A", "gsw")]
#exclude first induction by removing obs >= 547
inds <- inds[!(inds$plant == "IT82E-16_2" & inds$obs >= 547), ]

plot(A ~ Pci,
			data = AC[AC$plant == "IT82E-16_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = inds[inds$plant == "IT82E-16_9", ],
				col = 2
				)

lapply(levels(inds$geno), plot.Indonegeno, data = inds)

#close down all current graphics
graphics.off()

################################################################################
#update cp
cp <- rbind(AC, AQ, inds)
row.names(cp) <- c(1:nrow(cp))

lapply(levels(cp$geno), plot.Acionegeno, data = cp[cp$curve == "Aci", ])
lapply(levels(cp$geno), plot.AQonegeno, data = cp[cp$curve == "AQ", ])
lapply(levels(cp$geno), plot.Indonegeno, data = cp[cp$curve == "induction", ])

graphics.off()

