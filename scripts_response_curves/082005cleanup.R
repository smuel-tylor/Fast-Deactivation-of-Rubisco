#Produced for Nature Plants analysis to help de-clutter summarize.R scripts by making this data formatting and cleaning step operate within a single line
#updated 0820 for style
 
rm(list=ls())

library(lattice)

setwd("C:/Users/taylor53/Lancaster University/Carmo Silva, Elizabete - RIPE_Lancaster/Cowpea_GasExchange/Rdata")
#setwd("~/Downloads")

#####################################################
#####################################################
#load up and some basic housekeeping

cp <- read.csv("071819cowpeaGE.csv")
#rename genotypes and order appropriately for presentation with cultivated/uncultivated
geno.rename <- as.character(cp$geno)
geno.rename <- replace(geno.rename, geno.rename == "TVNu1948", "TVNu-1948")
geno.rename <- replace(geno.rename, geno.rename == "UIUC", "IT86D-1010")
geno.rename <- replace(geno.rename, geno.rename == "PI582537", "IT82E-16")
geno.rename <- replace(geno.rename, geno.rename == "PI582627", "KVu 379-P1")
cp$geno.old <- cp$geno
cp$geno <- factor(geno.rename,
									levels = c("Vadenantha", "TVNu-1948", "IT86D-1010", "IT82E-16", "KVu 379-P1", "Iron&Clay")
									)

#add a factor for cultivation status
cp$cult <- replace(as.character(cp$geno),
										cp$geno == "Vadenantha" | cp$geno == "TVNu-1948",
										"wild"
										)
cp$cult <- replace(cp$cult, cp$cult != "wild", "cultivar")
cp$cult <- factor(cp$cult)

#identifier for individuals
cp$plant <- factor(paste(cp$geno, cp$block, sep = "_"))
levels(cp$plant)

######################################################
######################################################
#preliminary data cleaning

x11()
xyplot(A ~ Ci | geno, data = cp[cp$curve == "Aci", ])
x11()
xyplot(A ~ Ci, groups = geno, auto.key = T, data = cp[cp$curve == "Aci", ])
#in analysis need to identify operating points - first value in each curve
#also need to clean out non-steady state values for 430 ppm (might be steady state for A,  but PhiPS2 often not equivalent to operating value

#plotting to evaluate errors in dataset and identify duplicated measurements
for (i in levels(cp$geno)){
	x11(w = 10, h = 6)
	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
	for ( j in as.numeric(levels(factor(cp$block))) ){
		dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
		if (nrow(dat) > 1){
			dat <- dat[order(dat$Pci), ]
			plot(A ~ Pci, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
			points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
			axis(side = 4, at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4), labels = c(0, .1, .2, .3, .4))
			mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
			} else {
					plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
					}
	}
}

#there are no clear issues here
#but need to take take care of the 'return to normal' points added to the end of curves
cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, "co2_at"]
#in an unsorted dataset it is usually points 8 and 17

#finter out these points and plot result
rownames(cp) <- c(1:nrow(cp))

for (i in levels(cp$geno)){
	x11(w = 10, h = 6)
	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
	for ( j in as.numeric(levels(factor(cp$block))) ){
	dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
	print(dat)
	if (nrow(dat) > 1){

		#filter out additional 430s
		cull <- rownames(dat)[c(8, 17)]
		cull <- cull[!is.na(cull)]
		for ( k in 1:length(cull) ){
			cp <- cp[rownames(cp)!= cull[k], ]
		 }
		rownames(cp) <- c(1:nrow(cp))
		#redo dat after filter
		dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
		#sort
		dat <- dat[order(dat$Pci), ]
		plot(A ~ Pci, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
		points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
		axis(side = 4,
					at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4),
					labels = c(0, .1, .2, .3, .4)
					)
		mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
		} else {
						plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
						}
		}
 }

#based on plots: check TVNu-1948-8
#(separate checks showed that I added a 1200 ppm point to this one,  which is why the final 430 wasn't trimmed (1200 was instead)
#I think it's fine to drop that data,  as it deviates from the protocol anyway
cp[cp$curve == "Aci" & cp$geno == "TVNu1948" & cp$block == 8, ]
cp <- cp[rownames(cp)!= 12943, ]

for (i in levels(cp$geno)){
	x11(w = 10, h = 6)
	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
	for ( j in as.numeric(levels(factor(cp$block))) ){
		dat <- cp[cp$curve == "Aci" & cp$geno == i & cp$block == j, ]
		print(dat)
		if (nrow(dat) > 1){
			dat <- dat[order(dat$Pci), ]
			plot(A ~ Pci, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
			points(110 * PhiPS2 ~ Pci, data = dat, pch = 21, bg = 2)
			axis(side = 4,
						at = c(0, 110 * 0.1, 110 * 0.2, 110 * 0.3, 110 * 0.4),
						labels = c(0, .1, .2, .3, .4)
						)
			mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
			} else {
							plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
							}
		}
	}
#GOOD!

for (i in levels(cp$geno)){
	x11(w = 10, h = 6)
	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
	for ( j in as.numeric(levels(factor(cp$block))) ){
		dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
		if (nrow(dat) > 1){
			dat <- dat[order(dat$Pci), ]
			plot(A ~ Qin, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
			points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
			axis(side = 4,
						at = c(0, 55 * 1 / 4, 55 * 1 / 2, 55 * 3 / 4, 55),
						labels = c(0, .25, .5, .75, 1)
						)
			mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
			} else {
							plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
							}
	}
}
#need to remove return to 1200 points from end of curves

rownames(cp) <- c(1:nrow(cp))

for (i in levels(cp$geno)){
	x11(w = 10, h = 6)
	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
	for ( j in as.numeric(levels(factor(cp$block))) ){
		dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
		print(dat)
		if (nrow(dat) > 1){

			#filter out additional 430s
			cull <- rownames(dat)[nrow(dat)]
			if (!is.na(cull)){
				for ( k in 1:length(cull) ){
					cp <- cp[rownames(cp) != cull, ]
					}
				rownames(cp) <- c(1:nrow(cp))
				#redo dat after filter
				dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
				#sort
				dat <- dat[order(dat$Qin), ]
				plot(A ~ Qin, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
				points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
				axis(side = 4,
							at = c(0, 55 * 1 / 4, 55 * 1 / 2, 55 * 3 / 4, 55),
							labels = c(0, .25, .5, .75, 1)
							)
				mtext(expression(Phi[PSII]), side = 4, line = 3, at = (60 / 2) - 5, las = 3)
				}
			} else {
							plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
							}
		}
	}


for (i in levels(cp$geno)){
	x11(w = 10, h = 6)
	par(mfrow = c(2, 4), mar = c(5, 4.5, 2.5, 4.5), las = 1)
	#par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 6), las = 1)
	for ( j in as.numeric(levels(factor(cp$block))) ){
		dat <- cp[cp$curve == "AQ" & cp$geno == i & cp$block == j, ]
		if (nrow(dat) > 1){
			dat <- dat[order(dat$Pci), ]
			plot(A ~ Qin, data = dat, main = paste(i, j, sep = " "), ylim = c(-5, 55))
			points(55 * PhiPS2 ~ Qin, data = dat, pch = 21, bg = 2)
			axis(side = 4,
						at = c(0, 55 * 1 / 4, 55 * 1 / 2, 55 * 3 / 4, 55),
						labels = c(0, .25, .5, .75, 1)
						)
			mtext(expression(Phi[PSII]), side = 4, line = 3, at = ( 60 / 2 ) - 5, las = 3)
						} else {
								plot(1, 1, type = "n", axes = F, bty = "n", main = j, xlab = "", ylab = "")
								}
		}
}

#look at IT86D-1010-9 and PI862537-2 inductions (cases where two inductions were recorded) - expect 2nd induction to be better
x11()
plot(gsw~induction.s,
			data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction", ], 
			ylim = c(0, 1),
			xlim = c(-1320, 1800)
			)
plot(A~induction.s,
			data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
plot(Pci~induction.s,
			data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)

cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "Aci", c("A", "gsw")]		

plot(A ~ Pci,
			data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "Aci", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction", ],
				col = 2
				)
#both responses overlay the A/ci response fairly well
#whilst the first response may better characterise induction without the modifications of gsw that result from Aci/AQ
#the second curve shows a more usual stomatal response during shade	
cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction", c("obs", "induction.s")]
#drop the first induction
#i.e.,  data with obs< = 363
cp <- cp[!(cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction" & cp$obs <= 363), ]

plot(A ~ Pci,
			data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "Aci", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp[cp$geno == "IT86D-1010" & cp$block == 9 & cp$curve == "induction", ],
				col = 2
				)

plot(gsw ~ induction.s,
			data = cp[cp$geno == "IT82E-16" & cp$block == 2 & cp$curve == "induction", ], 
			ylim = c(0, 1),
			xlim = c(-1320, 1800)
			)
#one of these is very noisy,  but in the other gsw does not recover

plot(A ~ induction.s,
			data = cp[cp$geno == "IT82E-16" & cp$block == 2 & cp$curve == "induction", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
plot(Pci ~ induction.s,
			data = cp[cp$geno == "IT82E-16" & cp$block == 2 & cp$curve == "induction", ], 
			ylim = c(0, 43),
			xlim = c(-1320, 1800)
			)
#it looks like the first attempt was better despite noisy ci
#this is because gsw recovers,  whereas
#second attempt does not show complete recovery
plot(A ~ Pci,
			data = cp[cp$geno == "IT82E-16" & cp$block == 9 & cp$curve == "Aci", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp[cp$geno == "IT82E-16" & cp$block == 9 & cp$curve == "induction", ],
				col = 2
				)

cp[cp$geno == "IT82E-16" & cp$block == 2 & cp$curve == "induction", c("obs", "induction.s", "A", "gsw")]
#exclude first induction by removing obs >= 547
cp <- cp[!(cp$geno == "IT82E-16" & cp$block == 2 & cp$curve == "induction" & cp$obs >= 547), ]

plot(A ~ Pci,
			data = cp[cp$geno == "IT82E-16" & cp$block == 9 & cp$curve == "Aci", ], 
			ylim = c(0, 43),
			xlim = c(0, 60)
			)
points(A ~ Pci,
				data = cp[cp$geno == "IT82E-16" & cp$block == 9 & cp$curve == "induction", ],
				col = 2
				)

#plot of all data was made at this stage, 
#after elimination of errors and removal of duplicate measurements
#this has not changed,  so code deleted here: original plot is "0818gasexchangeraw.pdf"

#close down all current graphics
graphics.off()
