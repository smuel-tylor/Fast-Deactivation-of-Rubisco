################################################################################
################################################################################
#For Supplement
#Redraw AC AJ analysis from nlmeVcmax.R using the modelled dataset
#This is parallel to 082005NaturePlantsnlmeVcmax_assessJlimitation.pdf
#but uses the subset of data identified above as of sufficient quality

library(here)
load(here("output/082005NaturePlantsnlmeVcmax.Rdata"))

#Re-subset the initial data.frame
#because the data used for modelling was a subset
subs_geno <- ilGKg.df$plant %in% unique(ilGKg.df.sub$plant)
ilGKg.df2 <- ilGKg.df[subs_geno, ]
ilGKg.df2$used <- ifelse(ilGKg.df2$induction.s > 1260 &
                           ilGKg.df2$induction.s <= 1500,
                         1,
                         0
)
ilGKg.df2$plant <- factor(gsub("TVNu-1948", "V. Sp. Savi", ilGKg.df2$plant))


#New helper functions that work with ilGKg.df.sub2
plot.AtAJ2 <- function(ind.AJ){
  plot(A ~ induction.s,
       data = ind.AJ,
       xlab = "time since end of shade (s)",
       ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1)),
       xlim = c(1200 - 120, 1200 + 900),
       ylim = c(0, 40),
       axes = FALSE,
       type = "n",
       main = ind.AJ$plant[1]
  )
  #line with transparency showing complete dataset
  lines(Aj.Pci ~ induction.s,
        data = ind.AJ,
        col = rgb(1, 0, 0, alpha = 0.3)
  )
  lines(A ~ induction.s,
        data = ind.AJ,
        col = rgb(0, 0, 0, alpha = 0.3)
  )
  #line without transparency showing modelled segment
  lines(Aj.Pci ~ induction.s,
        data = ind.AJ[ind.AJ$used == 1, ],
        col = rgb(1, 0, 0, alpha = 1),
        lwd = 2
  )
  lines(A ~ induction.s,
        data = ind.AJ[ind.AJ$used == 1, ],
        col = rgb(0, 0, 0, alpha = 1),
        lwd = 2
  )
  #reference line at 5 min
  lines(rep(1200 + 300, 2), c(-1e6, 1e6), lty = 2)
  
  axis(side = 1, at = seq(1200, 1200 + 900, 300), labels = seq(0, 900, 300))
  axis(side = 2, at = seq(0, 40, 20))
  
  box()
}

#Function looping the above.
#Here, input is AC.j modified induction data for multiple replicates
# in a dataframe, which is first converted to a list

#old version for reference
#plot.AtAJ.bygeno2 <- function(inds.AJ){
#  iAJ <- by
#  for (i in levels(inds$geno)){
#    par(mfrow = c(3, 2),
#        mar = c(4.5, 5, 2.5, 1.5),
#        mgp= c(3, 1, 0),
#        las = 1,
#        cex = 1
#    )
#    lapply(
#      inds.AJ[grep(i, names(inds.AJ), value = TRUE)],
#      plot.AtAJ
#    )
#  }
#}

plot.AtAJ.bygeno2 <- function(inds.AJ){
  for (i in levels(inds$geno)){
    par(mfrow = c(3, 2),
        mar = c(4.5, 5, 2.5, 1.5),
        mgp= c(2, 0.5, 0),
        las = 1,
        cex = 1,
        tcl = 0.2
    )
    iAJ <- inds.AJ[inds.AJ$geno == i, ]
    iAJ$plant <- factor(as.character(iAJ$plant))
    print(head(iAJ))
    by(
      iAJ,
      iAJ$plant,
      plot.AtAJ2
    )
  }
}

pdf(here("output/082005NaturePlantsFigureS8.pdf"),
    #paper = "a4",
    w = 15 / 2.58, h = 21.5 / 2.58
)

header.page(
  main = expression(
    atop("Does "*italic(A)*"(black line) intersect with",
         italic(A)[J] * " (red line),  predicted at steady state " *
           italic(J) * " & " * italic(c)[i] * " during induction"
    )
  )
)

plot.AtAJ.bygeno2(ilGKg.df2)

dev.off()
