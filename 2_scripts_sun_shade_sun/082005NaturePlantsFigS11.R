################################################################################
#Re-drawing the A/Q responses
# replicates nlmeLightResponses_nlme2.pdf
# but with TVNu renamed to V. sp.

library(here)
library(nlme)

load(here("output/082005NaturePlantsnlmeLightResponses.Rdata"))

################################################################################
#produce a plant-by-plant plot matching Vcmax and ActivationState
q.smth <- c(0:2000)
p.smth <- as.character(unique(AQ[ , "plant"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.AQ <- data.frame(
  Qin = rep(q.smth, length(p.smth)), 
  geno = factor(
    rep(g.smth, each = length(q.smth)),
    levels = levels(AQ$geno)
  ),
  plant = factor(
    rep(p.smth, each = length(q.smth)),
    levels = levels(AQ$plant)
  )
)

pred.AQ.nlme2 <- data.frame(
  A.geno = predict(AQ.nlme2, newdata = nd.AQ, level = 0),
  A.plant = predict(AQ.nlme2, newdata = nd.AQ, level = 1),
  plant = factor(nd.AQ$plant, levels = levels(nd.AQ$plant)),
  geno = factor(nd.AQ$geno, levels = levels(nd.AQ$geno)),
  Qin = nd.AQ$Qin 
)

fhds <- c("plant", "Qin", "A")
facs.AQ.nlme2 <- AQ[, fhds]
pred.AQ.nlme2 <- merge(pred.AQ.nlme2, facs.AQ.nlme2, all.x = TRUE)

#this is modified relative to nlmeLightResponses.R
levels(pred.AQ.nlme2$plant) <- gsub("TVNu-1948",
                                    "V. sp. Savi",
                                    levels(pred.AQ.nlme2$plant)
                                    )
levels(pred.AQ.nlme2$geno) <- gsub("TVNu-1948",
                                    "V. sp. Savi",
                                    levels(pred.AQ.nlme2$geno)
)

get_rep_list <-	function(geno, p.all){
  g.all <- p.all[p.all$geno == geno, ]
  levels(g.all$plant) <- replace(levels(g.all$plant),
                                 !levels(g.all$plant) %in% unique(g.all$plant),
                                 NA
  )
  by(g.all, g.all$plant, identity)
}

plot_rep <-	function(p.rep){
  ordp <- p.rep[order(p.rep$Qin), ]
  plot(1, 1,
       xlim = c(0, 2000),
       ylim = c(-3, 43),
       xlab = expression("incident PPFD " * (mu*mol~~m^-2~~s^-1)),
       ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1)),
       main = ordp$plant[1],
       type = "n",
       axes = FALSE)
  axis(side = 1, at = seq(0, 2000, 500), las = 1)
  axis(side = 2, at = seq(0, 40, 10), las = 1)
  points(ordp$A ~ ordp$Qin)
  lines(ordp$A.geno ~ ordp$Qin, lty = 1)
  lines(ordp$A.plant ~ ordp$Qin, lty = 2)
  box()
}

allrep_plot <- function(geno, p.all){
  p.reps <- get_rep_list(geno, p.all)
  par(mfrow = c(3, 2))
  lapply(p.reps, plot_rep)
}

pdf(here("output/082005NaturePlantsFigureS11.pdf"),
    w = 14.5 / 2.58, h = 22 / 2.58
)							

par(
  mar = c(4.5, 5.5, 3, 1),
  tcl = 0.4,
  oma = c(0, 0, 0, 0),
  cex = 1.2,
  cex.lab = 1.2,
  cex.axis = 1.2
)

glist <- levels(pred.AQ.nlme2$geno)
lapply(glist, allrep_plot, p.all = pred.AQ.nlme2)

dev.off()
