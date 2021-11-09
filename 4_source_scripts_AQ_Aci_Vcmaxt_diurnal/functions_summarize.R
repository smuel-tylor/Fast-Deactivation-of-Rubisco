#Functions written for/required for summarize.R

################################################################################
################################################################################

#General helper functions

#function that generates dataframes for output
cpdf <- function(vars){
  fcs <- c("plant", "geno", "block", "cult")
  lp <- levels(cp$plant)
  tdf <- data.frame(
    matrix(NA, length(lp), length(c(fcs, vars)))
  )
  names(tdf) <- c(fcs, vars)
  tdf$plant <- lp
  geno.strings <- unlist(strsplit(lp, split = "_"))
  tdf$geno <- factor(geno.strings[seq(1, length(geno.strings), 2)],
                     levels = levels(cp$geno)
  )
  tdf$block <- geno.strings[seq(2, length(geno.strings), 2)]
  tdf$cult <- replace(as.character(tdf$geno),
                      tdf$geno == "Vadenantha"|tdf$geno == "TVNu-1948",
                      "wild"
  )
  tdf$cult <- replace(tdf$cult, tdf$cult != "wild", "cultivar")
  tdf$cult <- factor(tdf$cult)
  tdf
}

#Insert blank page in pdf output containing a title
# (main is an expression or character vector)
header.page <- function(main){
  par(mfrow = c(1, 1), mar = c(10, 10, 10, 10))
  plot(1, 1,
       type = "n",
       bty = "n",
       axes = FALSE,
       xlab = "",
       ylab = "",
       main = main
  )
}

################################################################################
################################################################################

#functions for A/PPFD fitting

#implements AQfit, setting default values for the analysis in summarize.R
doAQfit <- function(AQ, 
                    params = c(phi = 0.05, Asat = 30, theta = 0.75, Rd = 1.9),
                    upp = NA,
                    low = NA,
                    Rd.fixed = FALSE
){
  if (!is.null(AQ)){
    AQ.cost.fits(AQ,
                 params = params,
                 upp = upp,
                 low = low,
                 Rd.fixed = Rd.fixed
    )
  }
}

#pars here is the naming wanted for the coefs, if this differs from 
mkAQsq <- function(AQfits, pars = NA){
  
  #function to extract summary parameters from ACi fits into a dataframe row
  # and rename as necessary
  pull.AQpars <- function(AQfit, pars){
    if (!is.null(AQfit)){
      out <- AQfit$coefs
      if (all(!is.na(pars))){
        names(out) <- pars
      }
      out
    } else { rep(NA, 4) }
  }
  
  AQmat <- t(sapply(AQfits, pull.AQpars, pars = pars))
  AQdf <- data.frame(apply(AQmat, 2, unlist))
  AQdf$plant <- row.names(AQdf)
  merge(cpdf(NULL), AQdf, all.x = TRUE)
}

#function to plot AQ fits as one page per genotype,
# two panels per plant showing whole fit and intercept region
plot.AQfit2 <- function(AQfit){
  if (!is.null(AQfit)){
    plot.AQfit(AQfit,
               xlim = c(0, 1600),
               ylim = c(-5, 55),
               PhiPSII = TRUE
    )
    title(main = AQfit$data$plant[1])
    plot.AQfit(AQfit,
               xlim = c(0, 150),
               ylim = c(-5, 20),
               PhiPSII = FALSE
    )
    title(main = AQfit$data$plant[1])
  }
}

plot.AQfit.bygeno <- function(AQfits){
  for (i in levels(AQ$geno)){
    par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
    lapply(AQfits[grep(i, names(AQfits), value = TRUE)], plot.AQfit2)
  }
}

################################################################################
################################################################################

#Equivalent functions for Aci fitting

#This function fits a curve using inputs and returns a fitted curve list object
#As above, the indexing and data selection are in a 'doACfits' function below
#List output means all elements of the fit are retained in these objects
#for subsequent plotting etc.
doACfit <- function(Aci, 
                    Rd.fixed = FALSE, 
                    Gamma.star = NA,
                    Gamma.star.fixed = FALSE, 
                    Kco = NA,
                    Kco.fixed = FALSE, 
                    gm = NA,
                    gm.fixed = FALSE
){
  
  #0121 added this for troubleshooting purposes
  print(as.character(Aci$plant[1]))
  A.op = Aci[1, "A"]
  Pci.op = Aci[1, "Pci"]
  Pca.op = Aci[1, "Pca"]
  
  Rd.guess <- Aci[1, "Rd.AQ"]
  if (is.na(Rd.guess)) {Rd.guess <- 2} #because Rd must be provided
  
  #generate the input
  #Gamma.star, Gamma.star.fixed, Kco, Kco.fixed, gm, gm.fixed, var.J, stoich
  input.AC <- Acifit.input(Aci, 
                           Rd = Rd.guess,
                           Rd.fixed, 
                           Gamma.star,
                           Gamma.star.fixed, 
                           Kco,
                           Kco.fixed, 
                           gm,
                           gm.fixed, 
                           stoich = "NADPH", 
                           var.J = FALSE
  )
  #add in lim states
  input.AC <- gen.lim.states(input.AC)
  #cost fits
  input.AC <- Aci.cost.fits(input.AC)
  #admissibility
  input.AC <- Acifit.admit(input.AC)
  #colimited points
  input.AC <- Acifit.colimit(input.AC)
  #choose model
  AdT <- input.AC$models[input.AC$models$admit == TRUE, ]
  best.AC = AdT[AdT$cost == min(AdT$cost, na.rm = TRUE), ]
  input.AC <- Acifit.evaluate(best.AC, input.AC, A.op, Pci.op, Pca.op)
  
  input.AC
}

#function to make a summary of key params from fits... complicated,
# but necessary to make a square df
# rewritten 0321 because there was re-ordering going on when this was produced
# not sure why mkAQsq ian't suffering from the same issues...
mkACsq <- function(ACfits){
  
  #parameters for Aci fits that are useful for summaries
  AC.pars = c("Vcmax", "J", "TPU", "Rd", "Gamma.star", "Kco", "gm",
              "A.40", "A.op", "Ls", "Pci.op", "Pca.op",
              "Gamma", "Ac.Aj", "Aj.Ap", "Ac.Ap",
              "Vcmax.25", "J.25", "TPU.25", "Rd.25", "gm.25"
  )
  
  #extract summary parameters from ACi fits into a dataframe row
  pull.ACpars <- function(input.AC, pars){
    core.pars <- c("Vcmax", "J", "TPU", "Rd", "Gamma.star", "Kco", "gm")
    out <- c(input.AC$chosen.mod[ , core.pars], 
             input.AC$F_Slimitation, 
             input.AC$ci.transitions, 
             input.AC$temperature.normalised
    )
    names(out) <- pars
    unlist(out)
  }
  
  ACdf <- cpdf(NULL) 
  ACdf[ , AC.pars] <- t(sapply(ACdf$plant,
                               function(.){
                                 if(!is.null(ACfits[[.]])) {
                                   pull.ACpars(ACfits[[.]], AC.pars)
                                 } else { rep(NA, length(AC.pars)) }
                               }
  )
  )
  ACdf

}

#function to plot Aci fits as one page per genotype,  two panels per plant
# showing whole fit and Gamma region
plot.ACfit2 <- function(input.AC){
  
  plot.Acifit(input.AC,
              xlim = c(0, 90),
              ylim = c(-5, 55),
              PhiPSII = TRUE,
              F_Slimitation = TRUE
  )
  title(main = input.AC$data$plant[1])
  plot.Acifit(input.AC,
              xlim = c(0, 20),
              ylim = c(-5, 20),
              PhiPSII = FALSE,
              F_Slimitation = FALSE
  )
  title(main = input.AC$data$plant[1])
  
}

plot.ACfit.bygeno <- function(AC.fits){
  for (i in levels(AC$geno)){
    par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
    lapply(AC.fits[grep(i, names(AC.fits), value = TRUE)], plot.ACfit2)
  }
}

################################################################################
################################################################################

#Astar and one-point Vcmax

#function that adds Astar and associated values to dataframe
do.Astarpredict <- function(inds, Aci){
  
  Aci.mod = Aci
  stoich = "NADPH"
  A = inds$A
  Pci = inds$Pci
  Pci.ss = Aci$Pci.op
  
  inds[, Astar.pars] <- one.point.Vcmax(Aci.mod, stoich, A, Pci, Pci.ss)
  inds
}

#plotting functions
plot.Astar2 <- function(inds, input.AC){
  
  plot.Acifit(input.AC,
              xlim = c(0, 90),
              ylim = c(-5, 55),
              phiPSII = TRUE,
              F_Slimitation = TRUE
  )
  points(A~Pci, data = inds, col = 3, pch = 21)
  title(main = input.AC$data$plant[1])
  
  plot(A~induction.s,
       data = inds, 
       xlim = c(-60, 900),
       ylim = c(0, 40), 
       xlab = expression(time~~since~~shade~~(s)), 
       ylab = expression(italic(A)~~(mu * mol~~m^-2~~s^-1)), 
       las = 1
  )
  points(Astar~induction.s, data = inds, col = 2, pch = 21)
  legend(900, 0,
         bty = "n",
         xjust = 1,
         yjust = 0,
         pch = rep(21, 2),
         col = c(1, 2),
         legend = expression(italic(A), italic(c)[i] *
                               " corrected,  " * italic(A) * " * "
                             )
  )
  title(main = input.AC$data$plant[1])
  
}

plot.Astar.bygeno <- function(inds.list, input.AC.list){
  for (i in levels(cp$geno)){
    par(mfrow = c(4, 2), mar = c(5, 6, 2.5, 1), las = 1)
    mapply(plot.Astar2,
           inds = inds.list[grep(i, names(inds.list), value = TRUE)],
           input.AC = input.AC.list[grep(i, names(input.AC.list), value = TRUE)]
    )
  }
}

#to make the above plot functions work,
# find the sets that have both an induction and an Aci response fitted
trim.inds <- function(inds.list, input.AC.list){
  if (length(input.AC.list) < length(inds.list)){
    inds.list <- inds.list[names(input.AC.list)]
  }
  
  if (length(input.AC.list) > length(inds.list)){
    input.AC.list <- input.AC.list[names(inds.list)]
  }
  
  inll <- !unlist(lapply(inds.list, is.null))
  anll <- !unlist(lapply(input.AC.list, is.null))
  
  inds.list[inll & anll]
  
}

trim.AC <- function(input.AC.list, inds.list){
  if (length(inds.list) < length(input.AC.list)){
    input.AC.list <- input.AC.list[names(inds.AC.list)]
  }
  
  if (length(inds.list) > length(input.AC.list)){
    inds.list <- inds.list[names(input.AC.list)]
  }
  
  inll <- !unlist(lapply(inds.list, is.null))
  anll <- !unlist(lapply(input.AC.list, is.null))
  
  input.AC.list[inll & anll]
}
