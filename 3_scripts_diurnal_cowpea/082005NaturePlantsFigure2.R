#Diagramming the diurnal model
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
                       alpha = 0.2 * 255
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

#load up diurnal models
#(includes light response model objects)
load(here("output/082005NaturePlantsdiurnalmodels.Rdata"))
objects()

#don't need to update geno names, as I'll just use one genotype to sketch the model
#... IT86D-1010, which has a particularly clear effect
##BEware if you try to adapt the below for TVNu! - see Figure1.R

#Nplants column width is ~8.9 cm = 3.44'
#inward ticks, sans serif
pdf(here("output/EGEetal_Fig2.pdf"), w = 3.5, h = 9)

#Activation State
par(
  mfrow = c(3, 1),
  mar = c(4, 6, 1, 2),
  oma = c(0, 0, 0, 0),
  mgp = c(3, 0.5, 0),
  tcl = 0.3,
  cex.axis = 1.2,
  cex.lab = 1.2,
  las = 1
)

#a - level 2 diurnal PPFD
plot(PPFD ~ Time.h, data = level2,
     type = "l",
     #lwd = 2,
     ylim = c(0, 2000),
     xlim = c(4.5, 19.5),
     ylab = expression(PPFD~~(mu*mol~~m^-2~~s^-1)),
     xlab = expression(time~~of~~day~~(hh:mm)),
     xaxt = "n"
)

axis(side = 1,
     at = c(6, 12, 18),
     labels = c("06:00", "12:00", "18:00"),
     lwd.ticks = 1,
     las = 1
)

mtext(expression(bold(a)),
      side = 2,
      at = 1.05 * 2000,
      line = 4.5,
      adj = 1
)

#b - detail of Aft and predict.A from midday period (11:45-13:30)
plot(predict ~ Time.h, data = dmods.Vc$"IT86D-1010",
     type = "n",
     ylim = c(0, 45),
     xlim = c(12.24, 14.24),
     ylab = expression(italic(A)~~(mu * mol~~m^-2~~s^-1)),
     xlab = expression(time~~of~~day~~(hh:mm)),
     xaxt = "n"
)

axis(side = 1,
     at = c(12.5, 13, 13.5, 14),
     labels = c("12:30", "13:00", "13:30", "14:00"),
     las = 1
)

mtext(expression(bold(b)), side = 2, at = 1.05 * 45, line = 4.5, adj = 1)

lines(predict.A + F.A ~ Time.h, data = dmods.Vc$"IT86D-1010",
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 1
)

lines(predict.A ~ Time.h, data = dmods.AS$"IT86D-1010",
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 5
)

lines(predict.A ~ Time.h, data = dmods.Vc$"IT86D-1010",
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 6
)


legend(12.24, 47,
       xjust = 0,
       yjust = 1,
       cex = 1.2,
       legend = expression(
         "steady-state PPFD response",
         italic(tau)["d,S"]~~italic(tau)["a,S"],
         italic(tau)["d,V"]~~italic(tau)["a,V"]
       ),
       lty = c(1,5,6),
       lwd = 1,
       col = o.cols["IT86D-1010"],
       ncol = 1,
       bty = "n"
)

#c - cumulative A using same data as c
#dmods objects that get rid of NA
dmods.Vc.nona <- lapply(dmods.Vc, function(.){ .[!is.na(.$Aft), ] })
dmods.AS.nona <- lapply(dmods.AS, function(.){ .[!is.na(.$Aft), ] })
dmods.VA.nona <- lapply(dmods.VA, function(.){ .[!is.na(.$Aft), ] })
dmods.AV.nona <- lapply(dmods.AV, function(.){ .[!is.na(.$Aft), ] })

ult <- dmods.Vc.nona$"IT86D-1010"$Time.h <= 19.5

plot(cumsum(Aft) / 1000 ~ Time.h,
     data = dmods.Vc.nona$"IT86D-1010"[ult, ],
     type = "n",
     ylim = c(0, 630),
     xlim = c(4.5, 19.5),
     ylab = expression(italic(A)~~(mmol~~m^-2)),
     xlab = expression(time~~of~~day~~(hh:mm)),
     xaxt = "n",
     yaxt = "n"
)

axis(side = 1,
     at = c(6, 12, 18),
     labels = c("06:00", "12:00", "18:00"),
     las = 1
)

axis(side = 2,
     at = c(0, 200, 400, 600),
     las = 1
)

mtext(expression(bold(c)), side = 2, at = 1.05 * 700, line = 4.5, adj = 1)

lines(cumsum(Aft) / 1000 ~ Time.h,
      data = dmods.Vc.nona$"IT86D-1010"[ult, ],
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 1
)

lines(cumsum(predict) / 1000 ~ Time.h,
      data = dmods.AS.nona$"IT86D-1010"[ult, ],
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 5
)

lines(cumsum(predict) / 1000 ~ Time.h,
      data = dmods.Vc.nona$"IT86D-1010"[ult, ],
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 6
)

lines(cumsum(predict) / 1000 ~ Time.h,
      data = dmods.VA.nona$"IT86D-1010"[ult, ],
      col = o.cols["IT86D-1010"],
      lwd = 1,
      lty = 2
)

lines(cumsum(predict) / 1000 ~ Time.h,
      data = dmods.AV.nona$"IT86D-1010"[ult, ],
      col = o.cols["IT86D-1010"],
      lwd = 2,
      lty = 3
)

legend(4.5, 660,
       xjust = 0,
       yjust = 1,
       cex = 1.2,
       legend = expression(
         "steady-state PPFD response",
         italic(tau)["d,S"]~~italic(tau)["a,S"],
         italic(tau)["d,V"]~~italic(tau)["a,V"],
         italic(tau)["d,V"]~~italic(tau)["a,S"],
         italic(tau)["d,S"]~~italic(tau)["a,V"]
       ),
       lty = c(1, 5, 6, 2, 3),
       lwd = c(rep(1, 4), 2),
       col = rep(o.cols["IT86D-1010"], 5),
       ncol = 1,
       bty = "n"
)

dev.off()

