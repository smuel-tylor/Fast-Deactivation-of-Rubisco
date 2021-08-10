################################################################################
#Plotting equivalent to nlme_4.pdf
#but with names matching those in the manuscript
library(here)
library(nlme)
load(here("output/082005NaturePlantsnlmeVcmax.Rdata"))

t.smth <- c(1200:1800)
p.smth <- as.character(unique(ilGKg.df.sub[ , "plant"]))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.Vc <- data.frame(induction.s = rep(t.smth, length(p.smth)), 
                    geno = rep(g.smth, each = length(t.smth)), 
                    plant = rep(p.smth, each = length(t.smth))
)
nd.Vc$plant <- factor(nd.Vc$plant, levels = levels(ilGKg.df.sub$plant))
nd.Vc$geno <- factor(nd.Vc$geno, levels = levels(ilGKg.df.sub$geno))

pred.ind.nlme4 <- data.frame(
  Vc.geno = predict(ind.nlme4, newdata = nd.Vc, level = 0), 
  Vc.plant = predict(ind.nlme4, newdata = nd.Vc, level = 1), 
  plant = nd.Vc$plant,
  geno = nd.Vc$geno, 
  induction.s = nd.Vc$induction.s
)

fhds <- c("plant", "induction.s", "Vcmax.t")
facs.ind.nlme4 <- ilGKg.df.sub[, fhds]
pred.ind.nlme4 <- merge(pred.ind.nlme4, facs.ind.nlme4, all.x = TRUE)

head(pred.ind.nlme4)
pred.ind.nlme4$plant <- factor(
  gsub("TVNu-1948", "V. sp. Savi", pred.ind.nlme4$plant)
)
pred.ind.nlme4$geno <- factor(
  gsub("TVNu-1948", "V. sp. Savi", pred.ind.nlme4$geno)
)

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
       ylim = c(0, 400),
       xlab = "time from end of shade (s)",
       ylab = expression(italic(V)["c,max"]~~(mu * mol~~m^-2~~s^-1)),
       main = ordp$plant[1],
       type = "n",
       axes = FALSE)
  axis(side = 1,
       at = seq(1200, 1800, 120),
       labels = seq(1200, 1800, 120) - 1200,
       las = 1
  )
  axis(side = 2,
       at = seq(0, 400, 100),
       las = 1
  )
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

pdf(here("output/082005NaturePlantsFigureS10.pdf"),
    w = 15 / 2.58, h = 22.5 / 2.58
)							

par(mar = c(4.5, 5.5, 3, 1),
    tcl = 0.4,
    oma = c(0, 0, 0, 0),
    cex = 1.2,
    cex.lab = 1.2,
    cex.axis = 1.2
)

glist <- levels(pred.ind.nlme4$geno)
lapply(glist, allrep_plot, p.all = pred.ind.nlme4)

dev.off()
