#Plotting predictions for Rubisco limitation based on gas exchange and
# Rubisco activation state
#0121 updated for here() & R 4.x, and to match new consistent geno/plant naming
# and corrections made to Vcmax analysis

library(here)

#colour schemes: indexed directly with species names for convenience
o.cols <- c(
  "V. adenantha" = rgb(0.5 * 255,
                       0.5 * 255,
                       0.5 * 255,
                       maxColorValue = 255
                       ),
  "V. sp. Savi." = rgb(0,
                       0,
                       0,
                       maxColorValue = 255
                       ),
  "IT82E-16" = rgb(135,
                   240,
                   115,
                   maxColorValue = 255
                   ),
  "IT86D-1010" = rgb(70,
                     190,
                     110,
                     maxColorValue = 255
                     )
)

t.cols <- c(
  "V. adenantha" = rgb(0.5 * 255,
                       0.5 * 255,
                       0.5 * 255,
                       maxColorValue = 255,
                       alpha = 0.2*255
                       ),
  "V. sp. Savi."  = rgb(0,
                        0,
                        0,
                        maxColorValue = 255,
                        alpha = 0.2 * 255),
  "IT82E-16" = rgb(135,
                   240,
                   115,
                   maxColorValue = 255,
                   alpha = 0.2 * 255
                   ),
  "IT86D-1010" = rgb(70,
                     190,
                     110,
                     maxColorValue = 255,
                     alpha = 0.2 * 255
                     )
)

#load up Vcmax time-series
#Noting that for this there is only one timeseries for Vcmax: no geno effects
load(here("output/082005NaturePlantsnlmeVcmax.Rdata"))

#Activation State
load(here("output/082005NaturePlantsnlmeActivationState.Rdata"))

#Update genotype names
TVNu.sub <- function(.){
  levels(.$geno) == "TVNu-1948"
}

levels(Vutl.noAdb$geno)[TVNu.sub(Vutl.noAdb)] <- "V. sp. Savi."
levels(pred.Vutl.nlme4$geno)[TVNu.sub(pred.Vutl.nlme4)] <- "V. sp. Savi."

levels(ilGKg.df.sub$geno)[TVNu.sub(ilGKg.df.sub)] <- "V. sp. Savi."
levels(in4$geno)[TVNu.sub(in4)] <- "V. sp. Savi."

#Nplants column width is ~8.9 cm = 3.44'
#inward ticks, sans serif

#Activation State
pdf(here("output/EGEetal_Fig1.pdf"), w = 3.5, h = 6.5)
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

plot(ActivationState ~ time.s,
     data = Vutl.noAdb,
     type = "n",
     ylim = c(0, 100),
     xlim = c(-5 * 60, 40 * 60),
     ylab = expression(italic(S) * " (%)"),
     xlab = expression(time~~from~~end~~of~~shade~~(min)),
     xaxt = "n"
)

axis(side = 1,
     at = c(0, 10, 20, 30, 40) * 60,
     labels = c(-20, -10, 0, 10, 20),
     cex.axis = 1,
     lwd.ticks = 1,
     las = 1
)

axis(side = 2,
     at = seq(0, 100, 20),
     lwd.ticks = 1,
     las = 1
)

mtext(expression(bold(a)),
      side = 2,
      at = 105,
      line = 4,
      adj = 1
)

legend(1200 + 10 * 60, 20,
       xjust = 0.5,
       yjust = 0.5,
       legend = c(
         as.expression(bquote(italic(.(levels(ilGKg.df.sub$geno)[1])))),
         bquote(italic(.(levels(ilGKg.df.sub$geno)[2]))),
         bquote(plain(.(levels(ilGKg.df.sub$geno)[3]))),
         bquote(plain(.(levels(ilGKg.df.sub$geno)[4])))
       ),
       pch = 21,
       pt.bg = t.cols,
       pt.lwd = 0,
       lwd = 2,
       col = o.cols,
       ncol = 1,
       bty = "n",
       cex = 0.7
)

lapply(
  levels(Vutl.noAdb$geno),
  function(.){
    points(ActivationState ~ time.s,
           data = Vutl.noAdb[Vutl.noAdb$geno == .,],
           pch = 21,
           col = NA,
           bg = t.cols[.],
           cex = 0.7
    )
  }
)

lapply(
  levels(Vutl.noAdb$geno),
  function(.){
    use <- pred.Vutl.nlme4[pred.Vutl.nlme4$geno == ., ]
    use <- use[c(1:2581), ]
    lines(AS.geno ~ time.s,
          data = use,
          col = o.cols[.],
          lwd = 2
    )
  }
)

#Vcmax
plot(Vcmax.t ~ induction.s,
     #0321 fixed here & below to match nlmeVcmax.R updates
     data = ilGKg.df.sub,
     type = "n",
     ylim = c(0, 370),
     xlim = c(20 * 60, 30 * 60),
     ylab = expression(italic(V)["c,max"]~~(mu * mol~~m^-2~~s^-1)),
     xlab = expression(time~~from~~end~~of~~shade~~(min)),
     xaxt = "n",
     yaxt = "n"
)

axis(side = 1,
     at = c(20, 22, 24, 26, 28, 30) * 60,
     labels = c(0, 2, 4, 6, 8, 10),
     cex.axis = 1,
     lwd.ticks = 1,
     las = 1
)

axis(side = 2,
     at = seq(0, 400, 100),
     lwd.ticks = 1,
     las = 1
)

mtext(expression(bold(b)),
      side = 2,
      at = 1.05 * 370,
      line = 4,
      adj = 1
)

lapply(
  levels(ilGKg.df.sub$geno),
  function(.){
    points(Vcmax.t ~ induction.s,
           data = ilGKg.df.sub[
             ilGKg.df.sub$geno == . &
               ilGKg.df.sub$induction.s > 1260,
             ],
           pch = 21,
           col = NA,
           bg = t.cols[.],
           cex = 0.7
    )
  }
)

lapply(
  levels(ilGKg.df.sub$geno),
  function(.){
    use <- in4[in4$geno == ., ]
    use <- use[c(1:121), ]
    print(use)
    lines(Vcmax.t ~ induction.s,
          data = use,
          col = o.cols[.],
          lwd = 2
    )
  }
)

dev.off()
