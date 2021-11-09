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

#load up Vcmax time-series
#Noting that for this there is only one timeseries for Vcmax: no geno effects
load(here("output/nlmeVcmax.Rdata"))

#Activation State
load(here("output/nlmeActivationState.Rdata"))

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
pdf(here("output/Figure1.pdf"), w = 7, h = 6.5)
#x11(w=3.5,h=6.5)

par(
  mfrow = c(2, 2),
  mar = c(4, 5.5, 1, 0.5),
  oma = c(0, 0, 0, 0),
  mgp = c(2, 0.5, 0),
  tcl = 0.3,
  cex.axis = 1.1,
  cex.lab = 1.1,
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
         bquote(italic(.(substr(levels(ilGKg.df.sub$geno)[2], 1, 6))) * " Savi."),
         bquote(plain(.(levels(ilGKg.df.sub$geno)[3]))),
         bquote(plain(.(levels(ilGKg.df.sub$geno)[4])))
       ),
       #pch = 21,
       #pt.bg = t.cols,
       #pt.lwd = 0,
       #lwd = 2,
       #col = o.cols,
       fill = o.cols,
       ncol = 1,
       bty = "n"#,
       #cex = 0.9
)

lapply(
  levels(Vutl.noAdb$geno),
  function(.){
    points(ActivationState ~ time.s,
           data = Vutl.noAdb[Vutl.noAdb$geno == .,],
           pch = 21,
           col = NA,
           bg = t.cols[.],
           cex = 0.8
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
     ylim = c(0, 340),
     xlim = c(20 * 60, 30 * 60),
     ylab = expression(italic(V)["c,max"]~~(mu * mol~~m^-2~~s^-1)),
     xlab = expression(time~~from~~end~~of~~shade~~(min)),
     xaxt = "n",
     yaxt = "n"
)

axis(side = 1,
     at = c(20, 22, 24, 26, 28, 30) * 60,
     labels = c(0, 2, 4, 6, 8, 10),
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
      at = 1.05 * 340,
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
           cex = 0.8
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

#Plotting Vi150 versus Vi850 and Vt150 versus Vt850

#subset to extract high and low PPFD data
#use the last three points at end of shade and end of sun
#because pre-shade is sometimes missing
#and this ensures n >= 2

Vutl.HL <- Vutl.noAdb[Vutl.noAdb$time.s >= 2040 |
                        Vutl.noAdb$time.s >= 725 &
                        Vutl.noAdb$time.s < 1200, ]

Vutl.HL$Vi_umol_m2s <- Vutl.HL$Vi_umol_min_mL *
  (((250 / 1000) / 0.55) * 1e4) / 60

Vutl.HL$Vt3_umol_m2s <- Vutl.HL$Vt3_umol_min_mL *
  (((250 / 1000) / 0.55) * 1e4) / 60

#where
# 250/1000 is the total extract volume in uL expressed as a fraction of a mL
# 0.55 cm^2 is the leaf disc area
# 1e4 multiplies up from cm^-2 to m^-2
# and the min^-1 rate is divided by 60 to obtain s^-1

Vutl.HL$PPFD <- ifelse(Vutl.HL$time.s >= 2040, 850, 150)

Vutl.H <- Vutl.HL[Vutl.HL$PPFD == 850, ]
Vutl.H.mn <- aggregate(Vutl.H[ , c("Vi_umol_m2s", "Vt3_umol_m2s")],
                       by = Vutl.H[, c("geno", "plant")],
                       mean, na.rm = TRUE
)
Vutl.H.sd <- aggregate(Vutl.H[ , c("Vi_umol_m2s", "Vt3_umol_m2s")],
                       by = Vutl.H[, c("geno", "plant")],
                       sd, na.rm = TRUE
)
Vutl.H.n <- aggregate(Vutl.H[ , c("Vi_umol_m2s", "Vt3_umol_m2s")],
                      by = Vutl.H[, c("geno", "plant")],
                      function(.){ length(.[!is.na(.)]) }
)
Vutl.H.summ <- merge(Vutl.H.mn, Vutl.H.sd,
                     by = c("geno", "plant"),
                     suffixes = c("_mn", "_sd")
)
Vutl.H.summ <- merge(Vutl.H.summ, Vutl.H.n,
                     by = c("geno", "plant")
)
names(Vutl.H.summ)[grep(".+m2s$", names(Vutl.H.summ))] <- paste0(
  names(Vutl.H.summ)[grep(".+m2s$", names(Vutl.H.summ))],
  "_n"
)
Vutl.L <- Vutl.HL[Vutl.HL$PPFD == 150, ]
Vutl.L.mn <- aggregate(Vutl.L[ , c("Vi_umol_m2s", "Vt3_umol_m2s")],
                       by = Vutl.L[, c("geno", "plant")],
                       mean, na.rm = TRUE
)
Vutl.L.sd <- aggregate(Vutl.L[ , c("Vi_umol_m2s", "Vt3_umol_m2s")],
                       by = Vutl.L[, c("geno", "plant")],
                       sd, na.rm = TRUE
)
Vutl.L.n <- aggregate(Vutl.L[ , c("Vi_umol_m2s", "Vt3_umol_m2s")],
                      by = Vutl.L[, c("geno", "plant")],
                      function(.){ length(.[!is.na(.)]) }
)
Vutl.L.summ <- merge(Vutl.L.mn, Vutl.L.sd,
                     by = c("geno", "plant"),
                     suffixes = c("_mn", "_sd")
)
Vutl.L.summ <- merge(Vutl.L.summ, Vutl.H.n,
                     by = c("geno", "plant")
)
names(Vutl.L.summ)[grep(".+m2s$", names(Vutl.L.summ))] <- paste0(
  names(Vutl.L.summ)[grep(".+m2s$", names(Vutl.L.summ))],
  "_n"
)

Vutl.HL.pairs <- merge(Vutl.H.summ, Vutl.L.summ,
                       by = c("geno", "plant"),
                       suffixes = c("_H", "_L")
) 

#check and update names
Vutl.HL.pairs$geno
#plant names also need updating for the following code
levels(Vutl.HL.pairs$plant) <- gsub("TVNu-1948",
                                    "V. sp. Savi.",
                                    levels(Vutl.HL.pairs$plant))
levels(Vutl.HL.pairs$plant)[
  !(levels(Vutl.HL.pairs$plant) %in% unique(Vutl.HL.pairs$plant))
] <- NA                                            

#add panels to figure
yupp <- 50

plot(Vi_umol_m2s_mn_L ~ Vi_umol_m2s_mn_H,
     data = Vutl.HL.pairs,
     type = "n",
     ylim = c(0, yupp),
     xlim = c(0, yupp),
     ylab = expression(italic(V)["i,150"]~~(mu*mol~~m^-2~~s^-1)),
     xlab = expression(italic(V)["i,850"]~~(mu*mol~~m^-2~~s^-1)),
     xaxt = "n",
     yaxt = "n"
)

axis(side = 1,
     at = seq(0, yupp, 10),
     lwd.ticks = 1,
     las = 1
)

axis(side = 2,
     at = seq(0, yupp, 10),
     lwd.ticks = 1,
     las = 1
)

abline(0, 1, lty = 2)

eqn <- lm(Vi_umol_m2s_mn_L ~ Vi_umol_m2s_mn_H - 1,
          data = Vutl.HL.pairs
)

abline(0, coef(eqn))

text(yupp, 6,
     bquote(y == .(round(coef(eqn), 3))*x),
     adj = c(1, 0.5)
)

mtext(expression(bold(c)),
      side = 2,
      at = 1.05 * yupp,
      line = 4,
      adj = 1
)

#legend(0, yupp + 2,
#       xjust = 0,
#       yjust = 1,
#       legend = c(
#         as.expression(bquote(italic(.(levels(Vutl.HL.pairs$geno)[1])))),
#         bquote(italic(.(levels(Vutl.HL.pairs$geno)[2]))),
#         bquote(plain(.(levels(Vutl.HL.pairs$geno)[3]))),
#         bquote(plain(.(levels(Vutl.HL.pairs$geno)[4])))
#       ),
#       pch = 19,
#       col = o.cols,
#       ncol = 1,
#       bty = "n",
#       cex = 0.8
#)

lapply(
  levels(Vutl.HL.pairs$plant),
  function(.){
    useit <- Vutl.HL.pairs[Vutl.HL.pairs$plant == ., ]
    lines(c(useit$Vi_umol_m2s_mn_H, useit$Vi_umol_m2s_mn_H),
          c(useit$Vi_umol_m2s_mn_L + useit$Vi_umol_m2s_sd_L,
            useit$Vi_umol_m2s_mn_L - useit$Vi_umol_m2s_sd_L
          ),
          col = o.cols[
            strsplit(.[1], "_")[[1]][1]
          ]
    )
    lines(c(useit$Vi_umol_m2s_mn_H + useit$Vi_umol_m2s_sd_H,
            useit$Vi_umol_m2s_mn_H - useit$Vi_umol_m2s_sd_H    
    ),
    c(useit$Vi_umol_m2s_mn_L,
      useit$Vi_umol_m2s_mn_L
    ),
    col = o.cols[
      strsplit(.[1], "_")[[1]][1]
    ]
    )
  }
)

lapply(
  levels(Vutl.HL.pairs$geno),
  function(.){
    useit <- Vutl.HL.pairs[Vutl.HL.pairs$geno == .,]
    points(Vi_umol_m2s_mn_L ~ Vi_umol_m2s_mn_H,
           data = useit,
           pch = 19,
           col = o.cols[.],
           cex = 0.8
    )
  }
)

#Vt
plot(Vt3_umol_m2s_mn_L ~ Vt3_umol_m2s_mn_H,
     data = Vutl.HL.pairs,
     type = "n",
     ylim = c(0, yupp),
     xlim = c(0, yupp),
     ylab = expression(italic(V)["t,150"]~~(mu*mol~~m^-2~~s^-1)),
     xlab = expression(italic(V)["t,850"]~~(mu*mol~~m^-2~~s^-1)),
     xaxt = "n",
     yaxt = "n"
)

axis(side = 1,
     at = seq(0, yupp, 10),
     lwd.ticks = 1,
     las = 1
)

axis(side = 2,
     at = seq(0, yupp, 10),
     lwd.ticks = 1,
     las = 1
)

mtext(expression(bold(d)),
      side = 2,
      at = 1.05 * yupp,
      line = 4,
      adj = 1
)

abline(0, 1, lty = 2)

eqnVt <- lm(Vt3_umol_m2s_mn_L ~ Vt3_umol_m2s_mn_H - 1,
            data = Vutl.HL.pairs
)

abline(0, coef(eqnVt))

text(yupp, 6,
     bquote(y == .(round(coef(eqnVt), 3))*x),
     adj = c(1, 0.5)
)


#legend(1200 + 10 * 60, 20,
#       xjust = 0.5,
#       yjust = 0.5,
#       legend = c(
#         as.expression(bquote(italic(.(levels(ilGKg.df.sub$geno)[1])))),
#         bquote(italic(.(levels(ilGKg.df.sub$geno)[2]))),
#         bquote(plain(.(levels(ilGKg.df.sub$geno)[3]))),
#         bquote(plain(.(levels(ilGKg.df.sub$geno)[4])))
#       ),
#       pch = 21,
#       pt.bg = t.cols,
#       pt.lwd = 0,
#       lwd = 2,
#       col = o.cols,
#       ncol = 1,
#       bty = "n",
#       cex = 0.7
#)

lapply(
  levels(Vutl.HL.pairs$plant),
  function(.){
    useit <- Vutl.HL.pairs[Vutl.HL.pairs$plant == ., ]
    lines(c(useit$Vt3_umol_m2s_mn_H, useit$Vt3_umol_m2s_mn_H),
          c(useit$Vt3_umol_m2s_mn_L + useit$Vt3_umol_m2s_sd_L,
            useit$Vt3_umol_m2s_mn_L - useit$Vt3_umol_m2s_sd_L
          ),
          col = o.cols[
            strsplit(.[1], "_")[[1]][1]
          ]
    )
    lines(c(useit$Vt3_umol_m2s_mn_H + useit$Vt3_umol_m2s_sd_H,
            useit$Vt3_umol_m2s_mn_H - useit$Vt3_umol_m2s_sd_H    
    ),
    c(useit$Vt3_umol_m2s_mn_L,
      useit$Vt3_umol_m2s_mn_L
    ),
    col = o.cols[
      strsplit(.[1], "_")[[1]][1]
    ]
    )
  }
)

lapply(
  levels(Vutl.HL.pairs$geno),
  function(.){
    useit <- Vutl.HL.pairs[Vutl.HL.pairs$geno == .,]
    points(Vt3_umol_m2s_mn_L ~ Vt3_umol_m2s_mn_H,
           data = useit,
           pch = 19,
           col = o.cols[.],
           cex = 0.8
    )
  }
)

dev.off()
