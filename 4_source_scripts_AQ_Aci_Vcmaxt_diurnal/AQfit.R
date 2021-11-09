#Fitting AQ responses using a an optim-based parallel with Acifit12
#only update from previous dated version is Theta -> theta
#as part of Git controlled changes for style
#- phiPSII -> PhiPSII
#- added ... in plot function
#- updated T and F to TRUE and FALSE, respectively

#Evaluate non-rectangular hyperbola model of AQ
#Q = PPFD, usually incident
#phi = initial slope (apparent quantum yield if using incident Q)
#Asat = asymptotic net CO2 assimilation at high Q
#theta = cruvature function determining how sharp the transition is
# between net CO2 assimilation at sub-saturating and saturating irradiances
#Rd = day assimilation: net CO2 assimilation when Q = 0 
AQ.form <- function(phi, Asat, theta, Rd, Q){
  (
    (phi * Q + Asat - sqrt((phi * Q + Asat)^2 - (4 * theta * phi * Q * Asat))) /
     (2 * theta)
   ) - Rd
	}

#AQ cost function
#i.e. find the minimum sum of squares
AQ.cost <- function(params, fixed = NA, data){
  
  if (is.na(fixed)) {
    phi = params["phi"]
    Asat = params["Asat"]
    theta = params["theta"]
    Rd = params["Rd"]
    } else {
      phi = params["phi"]
      Asat = params["Asat"]
      theta = params["theta"]
      Rd = fixed["Rd"]
      }
  
  Q = data[ , "Qin"]
  AQ.fitted = AQ.form(phi, Asat, theta, Rd, Q)
  AQ.resid = NA
  AQ.resid = data[ , "A"] - AQ.fitted
  
  sum(AQ.resid^2)
  
}

#Function using data and Aci.cost to fit the non-rectangular-hyperbola
#optionally takes fixed Rd
#params is a named vector with phi = , Amax = , theta = , Rd =
#upper and lower can optionally be set using named vectors too.
#If upper and lower are not set, they take default 'sensible' values
#doesn't currently report whether the model converged
AQ.cost.fits <- function(data, params, upp = NA, low = NA, Rd.fixed = FALSE){
  
  if (Rd.fixed == TRUE) {
    p = params[c("phi", "Asat", "theta")]
    f.p = params["Rd"]
    d.p = data[ , c("A", "Qin")]
    
    if (all(!is.na(upp))) {
      upp = upp
    } else {
        upp = c(phi = 0.2, Asat = 50, theta = 1)
    }
    
    if (all(!is.na(low))) {
      low = low
    } else {
        low = rep(1e-6, length(p))
        }
    
    } else {
      p = params[c("phi", "Asat", "theta", "Rd")]
      f.p = NA
      d.p = data[ , c("A","Qin")]
      
      if (all(!is.na(upp))) {
        upp = upp
      } else {
          upp = c(phi = 0.2, Asat = 50, theta = 1, Rd = 10)
      }
      
      if (all(!is.na(low))) {
        low = low
      } else {
          low = rep(1e-6, length(p))
      }
      }
  
  fit <- optim(p,
               fn = AQ.cost,
               fixed = f.p,
               data = d.p,
               method = "L-BFGS-B",
               upper = upp,
               lower = low,
               control = list(factr = 1e7, maxit = 500, trace = 1)
  )
  
  if (Rd.fixed == TRUE) {
    pars <- c(fit$par, f.p)
  } else {
      pars <- c(fit$par)
      }
  
  list(coefs = pars,
       start = p,
       upp = upp,
       low = low,
       Rd.fixed = Rd.fixed,
       data = data
  )
}

#plotting function with optional PhiPSII
#taking the output of AQ.cost.fits as input
plot.AQfit <- function(input, xlim = NULL, ylim = NULL, PhiPSII = FALSE, ...){
  
  Q <- c(0:2500)
  pAQ <- AQ.form(input$coefs["phi"],
                 input$coefs["Asat"],
                 input$coefs["theta"],
                 input$coefs["Rd"],
                 Q
                 )
	
	par(mar = c(5.5, 5.5, 3, 1), las = 1)
	
	if (is.null(ylim)) {
	  y.dat <- input$data[ , "A"]
	  ylim = range(y.dat[is.finite(y.dat)])
	  }
	
	if (is.null(xlim)) {
	  x.dat <- input$data[ , "Qin"]
	  xlim = range(x.dat[is.finite(x.dat)])
	  }
	
	if (PhiPSII == TRUE) {
	  par(mar = c(5.5, 5.5, 3, 5.5), las = 1)
	  }
	
	plot(A ~ Qin,
	     data = input$data,
	     xlim = xlim,
	     ylim = ylim,
	     type = "n",
	     xlab = expression(PPFD~~(mu*mol~~m^-2~~s^-1)),
	     ylab = expression(italic(A)~~(mu*mol~~m^-2~~s^-1))
	     )
	abline(0, 0)
	lines(pAQ ~ Q, lty = 2)
	points(A ~ Qin, data = input$data, pch = 21, bg = 1)
	
	if (PhiPSII == TRUE) {
	  phi.ex <- max(ylim) / (1.2 * max(input$data[ , "PhiPS2"], na.rm = TRUE))
	  points(phi.ex * PhiPS2 ~ Qin, data = input$data, pch = 24, bg = 0)
	  axis(side = 4,
	       at = seq(round(min(ylim), -1), round(max(ylim), -1), 10),
	       labels = round(
	         seq(round(min(ylim), -1), round(max(ylim), -1), 10) / phi.ex,
	         1
	         )
	       )
	  mtext(side = 4, expression(phi[PSII]), las = 3, line = 3, cex = 1.2)
	  legend(0.95 * diff(xlim) + min(xlim),
	         0.1 * diff(ylim) + min(ylim),
	         xjust = 1,
	         yjust = 0,
	         pch = c(21, 24),
	         pt.bg = c(1, 0),
	         lty = c(2, NA),
	         legend = expression(italic(A), italic(Phi)[PSII])
	  )
	}
}
