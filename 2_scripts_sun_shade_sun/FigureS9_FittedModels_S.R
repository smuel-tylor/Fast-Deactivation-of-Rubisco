################################################################################
#updating plot_nlme4.pdf with consistent naming
# specifically, renaming TVNu to V. sp

library(here)
library(nlme)

load(here("output/nlmeActivationState.Rdata"))

t.smth <- c(-240:2340)
s.0.smth <- ifelse(t.smth < 0,  1,  0)
s.d.smth <- ifelse(t.smth < 0 | t.smth >= 1200,  0,  1)
s.a.smth <- ifelse(t.smth >= 1200,  1,  0)
p.smth <- as.character(unique(Vutl.noAdb$plant))
g.smth <- sapply(strsplit(p.smth, "_"), function(.){ .[1] })
nd.AS <- data.frame(
  time.s = rep(t.smth, length(p.smth)), 
  s.0 = rep(s.0.smth, length(p.smth)), 
  s.d = rep(s.d.smth, length(p.smth)), 
  s.a = rep(s.a.smth, length(p.smth)), 
  geno = rep(g.smth, each = length(t.smth)), 
  plant = rep(p.smth, each = length(t.smth))
)
nd.AS$plant <- factor(nd.AS$plant, levels = levels(Vutl.noAdb$plant))
nd.AS$geno <- factor(nd.AS$geno, levels = levels(Vutl.noAdb$geno))

pred.Vutl.noAdb.nlme4 <- data.frame(
  AS.geno = predict(Vutl.nlme4, newdata = nd.AS, level = 0), 
  AS.plant = predict(Vutl.nlme4, newdata = nd.AS, level = 1), 
  plant = nd.AS$plant,
  geno = nd.AS$geno, 
  time.s = nd.AS$time.s 
)

fhds <- c("plant", "time.s", "ActivationState")
facs.Vutl.noAdb.nlme4 <- Vutl.noAdb[ , fhds]
pred.Vutl.noAdb.nlme4 <- merge(pred.Vutl.noAdb.nlme4,
                               facs.Vutl.noAdb.nlme4,
                               all.x = TRUE
                               )

#this is the key change relative to plot_nlme4.pdf
levels(pred.Vutl.noAdb.nlme4$plant) <- gsub("TVNu-1948",
                                            "V.sp. Savi",
                                            levels(pred.Vutl.noAdb.nlme4$plant)
)
levels(pred.Vutl.noAdb.nlme4$geno) <- gsub("TVNu-1948",
                                            "V.sp. Savi",
                                            levels(pred.Vutl.noAdb.nlme4$geno)
) 

get_rep_list <-	function(geno, p.all){
  g.all <- p.all[p.all$geno == geno, ]
  levels(g.all$plant) <- replace(
    levels(g.all$plant),
    !levels(g.all$plant) %in% unique(g.all$plant),
    NA
  )
  by(g.all, g.all$plant, identity)
}

plot_rep <-	function(p.rep){
  ordp <- p.rep[order(p.rep$time.s), ]
  plot(1, 1,
       xlim = c(-120, 2400),
       ylim = c(20, 100),
       xlab = "time from end of shade (s)", #"time from start of shade (s)",
       ylab = expression(italic(S)*" (%)"),
       main = ordp$plant[1],
       type = "n",
       axes = FALSE
  )
  axis(side = 1,
       at = seq(0, 2400, 600),
       labels = seq(0, 2400, 600) - 1200, #set axis origin to end of shade
       las = 1) 
  axis(side = 2,
       at = seq(0, 100, 20),
       las = 1)#, labels = rep("", length(seq(0, 100, 20))))
  points(ordp$ActivationState ~ ordp$time.s)
  lines(ordp$AS.geno ~ ordp$time.s, lty = 1)
  lines(ordp$AS.plant ~ ordp$time.s, lty = 2)
  box()
}

allrep_plot <- function(geno, p.all){
  p.reps <- get_rep_list(geno, p.all)
  par(mfrow = c(2, 2))
  lapply(p.reps, plot_rep)
}

pdf(here("output/FigureS9.pdf"),
    w = 6.5, h = 6.5
)

par(mar = c(4.5, 5.5, 3, 1), tcl = 0.4, oma = c(0, 0, 0, 0))

glist <- levels(pred.Vutl.noAdb.nlme4$geno)
lapply(glist, allrep_plot, p.all = pred.Vutl.noAdb.nlme4)

dev.off()

#looks OK,  and appears to justify the random effects being included
#tau.d and the activation state differ among genotypes,  tau.a is consistent
