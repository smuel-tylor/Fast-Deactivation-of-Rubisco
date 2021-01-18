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

#for convenience plotting 
library(lattice)
library(here)

################################################################################
#load up and basic housekeeping

cp <- read.csv(here("data/RIPE_20210118_cowpeaGE.csv"),
               stringsAsFactors = FALSE
               )

#rename genotypes and order appropriately based on cultivated/uncultivated
#next hashed lines written for R 3.x and now redundant 
#geno.rename <- as.character(cp$geno)
#geno.rename <- replace(geno.rename, geno.rename == "TVNu1948", "TVNu-1948")
#geno.rename <- replace(geno.rename, geno.rename == "UIUC", "IT86D-1010")
#geno.rename <- replace(geno.rename, geno.rename == "PI582537", "IT82E-16")
#geno.rename <- replace(geno.rename, geno.rename == "PI582627", "KVu 379-P1")

#these redundant with the new input file
#copy original names just in case 
#cp$geno.old <- cp$geno
#make names consistent/best
#cp$geno <- replace(cp$geno, cp$geno == "TVNu1948", "TVNu-1948")
#cp$geno <- replace(cp$geno, cp$geno == "UIUC", "IT86D-1010")
#cp$geno <- replace(cp$geno, cp$geno == "PI582537", "IT82E-16")
#cp$geno <- replace(cp$geno, cp$geno == "PI582627", "KVu 379-P1")

#add cultivation status
cp$cult <- replace(cp$geno,
                   cp$geno == "V. adenantha" | cp$geno == "TVNu-1948",
                   "wild"
                   )
cp$cult <- replace(cp$cult, cp$cult != "wild", "cultivar")

#make these into factors
g.levs <- c("V. adenantha", "TVNu-1948", "IT86D-1010", "IT82E-16")
#           , "KVu 379-P1", "Iron&Clay"

c.levs <- c("cultivar", "wild")

cp$geno <- factor(cp$geno, levels = g.levs)
cp$cult <- factor(cp$cult, levels = c.levs)

cp$block <- factor(cp$block)

#identifier for individuals
cp$plant <- factor(paste(cp$geno, cp$block, sep = "_"))
levels(cp$plant)

################################################################################
#preliminary data cleaning

x11()
xyplot(A ~ Ci | geno, data = cp[cp$curve == "Aci", ])
x11()
xyplot(A ~ Ci, groups = geno, auto.key = T, data = cp[cp$curve == "Aci", ])
#in analysis need to identify operating points - first value in each curve
#also need to clean out non-steady state values for 430 ppm (might be steady state for A,  but PhiPS2 often not equivalent to operating value

#plotting to evaluate errors in dataset and identify duplicated measurements

#this hashed out code replaced by below
##for (i in levels(cp$geno)){
##	x11(w = 10, h = 6)
##	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
##	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
##	for ( j in as.numeric(levels(factor(cp$block))) ){
##		dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
##		if (nrow(dat) > 1){
##			dat <- dat[order(dat$Pci), ]
##			plot(A ~ Pci, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
##			points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
##			axis(side = 4, at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4), labels = c(0, .1, .2, .3, .4))
##      mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
##			} else {
##					plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
##					}
##	}
##}

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

#plotAciifis(cp.Aci[cp.Aci$plant == levels(cp.Aci$plant)[1], ])

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

cp.Aci <- cp[cp$curve == "Aci", ]

lapply(levels(cp.Aci$geno), plot.Acionegeno, data = cp.Aci)

#there are no clear issues here

#but will remove the redundant 'return to normal' points
#because these result in unequal weighting of the dataset
#around the operating point, i.e., they represent pseudoreps

#an example
cp.Aci[cp.Aci$plant == "IT82E-16_1", "co2_at"]
#in an unsorted dataset it is usually points 8 and 17

#filter out these points and plot result
cp.Aci <- by(cp.Aci, cp.Aci$plant, identity)
cp.Aci <- lapply(cp.Aci, function(.){ .[c(1:7,9:16), ] })
cp.Aci <- do.call(rbind, cp.Aci)

lapply(levels(cp.Aci$geno), plot.Acionegeno, data = cp.Aci)

#below was replaced by the above
#for (i in levels(cp$geno)){
#	x11(w = 10, h = 6)
#	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
#
#	
#		for ( j in as.numeric(levels(factor(cp$block))) ){
#	dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
#	print(dat)
#	if (nrow(dat) > 1){
#
#		#filter out additional 430s
#		cull <- rownames(dat)[c(8, 17)]
#		cull <- cull[!is.na(cull)]
#		for ( k in 1:length(cull) ){
#			cp <- cp[rownames(cp)!= cull[k], ]
#		 }
#		rownames(cp) <- c(1:nrow(cp))
#		#redo dat after filter
#		dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
#		#sort
#		dat <- dat[order(dat$Pci), ]
#		plot(A ~ Pci, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
#		points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
#		axis(side = 4,
#					at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4),
#					labels = c(0, .1, .2, .3, .4)
#					)
#		mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
#		} else {
#						plot(1, 1,
#                  type = "n", axes = F,
#                  bty = "n", main = j, xlab = "", ylab = "")
#						}
#		}
# }

#the below doesn't apply because of the way the filtering now works
#based on plots: check TVNu-1948-8
#(separate checks showed that I added a 1200 ppm point to this one,
# which is why the final 430 wasn't trimmed (1200 was instead))
#I think it's fine to drop that data,  as it deviates from the protocol anyway
#cp[cp$geno == "TVNu-1948" & cp$block == 8, c(1:15)]
#cp <- cp[rownames(cp)!= 12943, ]

#for (i in levels(cp$geno)){
#	x11(w = 10, h = 6)
#	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
#	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
#	for ( j in as.numeric(levels(factor(cp$block))) ){
#		dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
#		print(dat)
#		if (nrow(dat) > 1){
#			dat <- dat[order(dat$Pci), ]
#			plot(A ~ Pci, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
#			points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
#			axis(side = 4,
#						at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4),
#						labels = c(0, .1, .2, .3, .4)
#						)
#			mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
#			} else {
#							plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
#							}
#		}
#	}
#GOOD!

#The below replaced as above
#for (i in levels(cp$geno)){
#	x11(w = 10, h = 6)
#	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
#	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
#	for ( j in as.numeric(levels(factor(cp$block))) ){
#		dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
#		if (nrow(dat) > 1){
#			dat <- dat[order(dat$Pci), ]
#			plot(A ~ Qin, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
#			points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
#			axis(side = 4,
#						at = c(0, 55 * 1 / 4, 55 * 1 / 2, 55 * 3 / 4, 55),
#						labels = c(0, .25, .5, .75, 1)
#						)
#			mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
#			} else {
#							plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
#							}
#	}
#}

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

#plotAQifis(cp.AQ[cp.AQ$plant == levels(cp$plant)[1], ])

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

cp.AQ <- cp[cp$curve == "AQ", ]

lapply(levels(cp.AQ$geno), plot.AQonegeno, data = cp.AQ)

#Remove 'return to 1500' points from end of curves

#an example
cp.AQ[cp.AQ$plant == "IT82E-16_1", "Qin"]
#in an unsorted dataset it is usually points 8 and 17

cp.AQ <- by(cp.AQ, cp.AQ$plant, identity)
cp.AQ <- lapply(cp.AQ,
                function(.){
                  if (length(.[, 1]) > 0) {
                    .[c(1:(length(.[, 1]) - 1)), ]
                  }
                  }
                )
cp.AQ <- do.call(rbind, cp.AQ)

lapply(levels(cp.AQ$geno), plot.AQonegeno, data = cp.AQ)

#rownames(cp) <- c(1:nrow(cp))

#for (i in levels(cp$geno)){
#	x11(w = 10, h = 6)
#	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
#	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
#	for ( j in as.numeric(levels(factor(cp$block))) ){
#		dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
#		print(dat)
#		if (nrow(dat) > 1){
#
#			#filter out additional 430s
#			cull <- rownames(dat)[nrow(dat)]
#			if (!is.na(cull)){
#				for ( k in 1:length(cull) ){
#					cp <- cp[rownames(cp) != cull, ]
#					}
#				rownames(cp) <- c(1:nrow(cp))
#				#redo dat after filter
#				dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
#				#sort
#				dat <- dat[order(dat$Qin), ]
#				plot(A ~ Qin, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
#				points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
#				axis(side = 4,
#							at = c(0, 55 * 1 / 4, 55 * 1 / 2, 55 * 3 / 4, 55),
#							labels = c(0, .25, .5, .75, 1)
#							)
#				mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
#				}
#			} else {
#							plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
#							}
#		}
#	}


#for (i in levels(cp$geno)){
#	x11(w = 10, h = 6)
#	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
#	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
#	for ( j in as.numeric(levels(factor(cp$block))) ){
#		dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
#		if (nrow(dat) > 1){
#			dat <- dat[order(dat$Pci), ]
#			plot(A ~ Qin, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
#			points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
#			axis(side = 4,
#						at = c(0, 55 * 1 / 4, 55 * 1 / 2, 55 * 3 / 4, 55),
#						labels = c(0, .25, .5, .75, 1)
#						)
#			mtext(expression(Phi[PSII]), side = 4, line = 3, at = ( 60 / 2 ) - 5,
#           las = 3
#           )
#						} else {
#								plot(1, 1, type = "n", axes = F,
#               bty = "n", main = j, xlab = "", ylab = "")
#								}
#		}
#}

#look at IT86D-1010-9 and PI862537-2 inductions
# (cases where two inductions were recorded)
# - expect 2nd induction to be better
cp.ind <- cp[cp$curve == "induction", ]

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

#plotIndifis(cp.ind[cp.ind$plant == levels(cp.ind$plant)[1], ])

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

lapply(levels(cp.ind$geno), plot.Indonegeno, data = cp.ind)

x11()
plot(gsw~induction.s,
			data = cp.ind[cp.ind$plant == "IT86D-1010_9", ],
			ylim = c(0, 1),
			xlim = c(-1320, 1800)
			)
plot(A~induction.s,
			data = cp.ind[cp.ind$plant == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
plot(Pci~induction.s,
			data = cp.ind[cp.ind$geno == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)

cp[cp$plant == "IT86D-1010_9" & cp$curve == "Aci", c("A", "gsw")]		

plot(A ~ Pci,
			data = cp.Aci[cp.Aci$plant == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp.ind[cp.ind$plant == "IT86D-1010_9", ],
				col = 2
				)
#Both responses overlay the A/ci response fairly well.
#The first response may better characterise induction without
# the modifications of gsw that result from Aci/AQ
#the second curve shows a more usual stomatal response during shade	
cp.ind[cp.ind$plant == "IT86D-1010_9" & cp.ind$curve == "induction",
   c("obs", "induction.s")]
#drop the first induction
#i.e.,data with obs <= 363
cp.ind <- cp.ind[!(cp.ind$plant == "IT86D-1010_9" & cp.ind$obs <= 363), ]

plot(A ~ Pci,
			data = cp.Aci[cp.Aci$plant == "IT86D-1010_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp.ind[cp.ind$plant == "IT86D-1010_9", ],
				col = 2
				)

plot(gsw ~ induction.s,
			data = cp.ind[cp.ind$plant == "IT82E-16_2", ], 
			ylim = c(0, 1),
			xlim = c(-1320, 1800)
			)
#one of these is very noisy,  but in the other gsw does not recover

plot(A ~ induction.s,
			data = cp.ind[cp.ind$plant == "IT82E-16_2", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
plot(Pci ~ induction.s,
			data = cp.ind[cp.ind$plant == "IT82E-16_2", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
#it looks like the first attempt was better despite noisy ci
#this is because gsw recovers,  whereas
#second attempt does not show complete recovery
plot(A ~ Pci,
			data = cp.Aci[cp$plant == "IT82E-16_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp.ind[cp.ind$plant == "IT82E-16_9", ],
				col = 2
				)

cp.ind[cp.ind$plant == "IT82E-16_2", c("obs", "induction.s", "A", "gsw")]
#exclude first induction by removing obs >= 547
cp.ind <- cp.ind[!(cp.ind$plant == "IT82E-16_2" & cp.ind$obs >= 547), ]

plot(A ~ Pci,
			data = cp.Aci[cp.Aci$plant == "IT82E-16_9", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp.ind[cp.ind$plant == "IT82E-16_9", ],
				col = 2
				)

lapply(levels(cp.ind$geno), plot.Indonegeno, data = cp.ind)

#close down all current graphics
graphics.off()

cp <- rbind(cp.Aci, cp.AQ, cp.ind)
row.names(cp) <- c(1:nrow(cp))

lapply(levels(cp$geno), plot.Acionegeno, data = cp[cp$curve == "Aci", ])
lapply(levels(cp$geno), plot.AQonegeno, data = cp[cp$curve == "AQ", ])
lapply(levels(cp$geno), plot.Indonegeno, data = cp[cp$curve == "induction", ])

graphics.off()
