#A new version of AQ fitting using a similar optim-based approach to the Acifit routines I have developed
#but obviously without the need for the limitation state matrices
#only update from prev version is Theta -> theta

#evaluates AQ
AQ.form<-function(phi,Asat,theta,Rd,Q){
	((phi*Q+Asat-sqrt((phi*Q+Asat)^2-(4*theta*phi*Q*Asat)))/(2*theta))-Rd
	}

#AQ cost function
AQ.cost<-function(params,fixed=NA,data){
	
	if (is.na(fixed)) { 	
							phi=params["phi"]
							Asat=params["Asat"]
							theta=params["theta"]
							Rd=params["Rd"] 
						} else {
								phi=params["phi"]
								Asat=params["Asat"]
								theta=params["theta"]
								Rd=fixed["Rd"] 
								}

	Q=data[,"Qin"]
	AQ.fitted=AQ.form(phi,Asat,theta,Rd,Q)
	AQ.resid=NA
	AQ.resid=data[,"A"]-AQ.fitted
	
	sum(AQ.resid^2)

}

#function that uses data and Aci.cost to fit model
#optionally takes fixed Rd
#params is named vector with phi,Amax,theta,Rd
#upper and lower can optionally be set using named vectors too
#ought to build in ability to see if model converged
AQ.cost.fits<-function(data,params,upp=NA,low=NA,Rd.fixed=F){

	if (Rd.fixed==T){ 	p=params[c("phi","Asat","theta")]
						f.p=params["Rd"]
						d.p=data[,c("A","Qin")]
						if (all(!is.na(upp))) { upp=upp } else { upp=c(phi=0.2,Asat=50,theta=1) }
						if (all(!is.na(low))) { low=low } else { low=rep(1e-6,length(p)) }
						} else {	p=params[c("phi","Asat","theta","Rd")]
									f.p=NA
									d.p=data[,c("A","Qin")]
									if (all(!is.na(upp))) { upp=upp } else { upp=c(phi=0.2,Asat=50,theta=1,Rd=10) }
									if (all(!is.na(low))) { low=low } else { low=rep(1e-6,length(p)) }
									}
	
	fit<-optim(p,
			fn=AQ.cost,
			fixed=f.p,data=d.p,
			method="L-BFGS-B",
			upper=upp,lower=low,
			control=list(factr=1e7,maxit=500,trace=1))

	if (Rd.fixed==T) { pars<-c(fit$par,f.p) } else { pars<-c(fit$par) }
	
	list(coefs=pars,start=p,upp=upp,low=low,Rd.fixed=Rd.fixed,data=data)
	
	}
	
plot.AQfit<-function(input,xlim=NULL,ylim=NULL,phiPSII=F){
	
	Q<-c(0:2500)
	pAQ<-AQ.form(input$coefs["phi"],input$coefs["Asat"],input$coefs["theta"],input$coefs["Rd"],Q)
	
	par(mar=c(5.5,5.5,3,1),las=1)
	
	if (is.null(ylim)) { y.dat<-input$data[,"A"]
										ylim=range(y.dat[is.finite(y.dat)]) }

	if (is.null(xlim)) { x.dat<-input$data[,"Qin"]
										xlim=range(x.dat[is.finite(x.dat)]) }
	
	if (phiPSII==T) {
				par(mar=c(5.5,5.5,3,5.5),las=1)
				}

	plot(A~Qin,data=input$data,xlim=xlim,ylim=ylim,type="n",
		xlab=expression(PPFD~~(mu*mol~~m^-2~~s^-1)),
		ylab=expression(italic(A)~~(mu*mol~~m^-2~~s^-1))
		)
	abline(0,0)
	lines(pAQ~Q,lty=2)
	points(A~Qin,data=input$data,pch=21,bg=1)

	if (phiPSII==T) {
				
					phi.ex<-max(ylim)/(1.2*max(input$data[,"PhiPS2"],na.rm=T))
					points(phi.ex*PhiPS2~Qin,data=input$data,pch=24,bg=0)
					axis(side=4,at=seq(round(min(ylim),-1),round(max(ylim),-1),10),labels=round(seq(round(min(ylim),-1),round(max(ylim),-1),10)/phi.ex,1))
					mtext(side=4,expression(phi[PSII]),las=3,line=3,cex=1.2)
					
					legend(0.95*diff(xlim)+min(xlim),0.1*diff(ylim)+min(ylim),xjust=1,yjust=0,
							pch=c(21,24),pt.bg=c(1,0),lty=c(2,NA),
							legend=expression(italic(A),italic(phi)[PSII])
							)
					}
	}

