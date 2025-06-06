#rm(list=ls(all=TRUE))

## ----setupBayes
restore<-TRUE
#it takes some time to run the analysis and the results differ slightly run-to-run
#use a previously created analysis osdibf.RData or start afresh by setting restore<-FALSE
session<-'osidbf.RData'
if(restore==TRUE){
    load(session)
} else {#create session
  
library(polspline)
library(BayesFactor)
library(R2jags)
library(BEST)
setDirs<-function(){#{{{
    cwd<-getwd()
    rstr<-'c:\\users\\admin\\'
    #rstr<-'c:\\users\\spg\\'
    setwd(paste0(rstr,'documents\\latex\\work\\osid'))
    source('..\\..\\..\\rwork\\functions\\generic.r')#assumes we are currently in latex\work\XXX
    library(Hmisc) 
    return(cwd)
}#}}}
cwd<-setDirs()
mydata<-read.table('../../../RWork/E175/itest.txt',header=FALSE,sep="\t",comment.char="")
colnames(mydata)<-c('sno','xrsp','igrp')


dx <- mydata$xrsp[mydata$igrp==0]
dy <- mydata$xrsp[mydata$igrp==1]
vt<-var.test(dx,dy)
varTest<-paste0('F(',vt$parameter[1],',',vt$parameter[2],')=',fn$fndp(vt$statistic,2),fn$formatp(vt$p.value))

n1 <- length(dx)
n2 <- length(dy)

# Rescale
#pooled sd 
nx<-length(dx)-1
ny<-length(dy)-1
psd<-((nx*var(dx)+ny*var(dy))/(nx+ny))^0.5
sy <- dy - mean(dx)#scaled x and y for model1
sy <- sy/psd
sx <- (dx-mean(dx))/psd; 
meanx<-mean(dx)

model1Data <- list('sx', 'sy', 'n1', 'n2') # to be passed on to JAGS
model2Data <- list('dx','dy','n1','n2','psd','meanx')
model3Data <- list('dx','dy','n1','n2','meanx')
	      

model1Inits <- list(
  list(delta = rnorm(1,0,2), mu = rnorm(1,0,2), sigmatmp = runif(1,0,5)),
  list(delta = rnorm(1,0,2), mu = rnorm(1,0,2), sigmatmp = runif(1,0,5)),
  list(delta = rnorm(1,0,2), mu = rnorm(1,0,2), sigmatmp = runif(1,0,5)))

model2Inits<-list(
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigma = runif(1,0,5)),
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigma = runif(1,0,5)),
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigma = runif(1,0,5)))

model3Inits<-list(
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigmatmp = runif(1,0,5)),
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigmatmp = runif(1,0,5)),
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigmatmp = runif(1,0,5)))

model4Inits<-list(
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigma = runif(1,0,5),vMinus=runif(1,3,7)),
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigma = runif(1,0,5),vMinus=runif(1,3,7)),
  list(delta = rnorm(1,0,2), mux = rnorm(1,0,2), sigma = runif(1,0,5),vMinus=runif(1,3,7)))

model1Parameters <- c('delta','mux','muy','sigma')
model2Parameters <- c('delta','mux','sigma') 
model4Parameters <- c('delta','mux','vMinus','sigma')

modelFile<-'temp.txt'

getModel<-function(m){#{{{
    if(m==1){#{{{
	result<-'model{ 
	      # Data
	      for (i in 1:n1){
		sx[i] ~ dnorm(mux,lambda)
	      }
	      for (j in 1:n2){
		sy[j] ~ dnorm(muy,lambda)
	      }
	      # Means and precision
	      alpha <- delta*sigma
	      mux <- mu	
	      muy <- mu+alpha
	      lambda <- pow(sigma,-2)
	      # delta, mu, and sigma Come From Cauchy Distributions
	      lambdadelta ~ dchisqr(1)
	      delta ~ dnorm(0,lambdadelta)
	      lambdamu ~ dchisqr(1)
	      mu ~ dnorm(0,lambdamu)
	      lambdasigma ~ dchisqr(1)
	      sigmatmp ~ dnorm(0,lambdasigma)
	      sigma <- abs(sigmatmp)#positive values only
	    }'
	}#}}}
     if(m==2){#{{{
	result<-'model{
	      for (i in 1:n1){
		dx[i] ~ dnorm(mux,lambda)
	      }
	      for (j in 1:n2){
		dy[j] ~ dnorm(mux+delta,lambda)
	      }
	      sigma~dunif(psd/100,psd*100)
	      lambda<-pow(sigma,-2)
	      mux~dnorm(meanx,1/(100*psd)^2)
	      delta~dnorm(0,1/(100*2^0.5*psd)^2)#assuming a common variance the variance of the difference is 2x common var -> 2^0.5 common sd
	}'
	}#}}}
     if(m==3){#{{{
	result<-'model{
	      for (i in 1:n1){
		dx[i] ~ dnorm(mux,lambda)
	      }
	      for (j in 1:n2){
		dy[j] ~ dnorm(mux+delta,lambda)
	      }
	      lambda<-pow(sigma,-2)
	      lambdasigma ~ dchisqr(1)
	      sigmatmp ~ dnorm(0,lambdasigma)
	      sigma <- abs(sigmatmp)#positive values only
	      lambdamux~dchisqr(1)
	      mux~dnorm(meanx,lambdamux)
	      lambdadelta~dchisqr(1)
	      delta~dnorm(0,lambdadelta)
	}'
	}#}}}
     if(m==4){#{{{
	result<-'model{
	      for (i in 1:n1){
		dx[i] ~ dt(mux,lambda,v)
	      }
	      for (j in 1:n2){
		dy[j] ~ dt(mux+delta,lambda,v)
	      }
	      sigma~dunif(psd/100,psd*100)
	      lambda<-pow(sigma,-2)
	      mux~dnorm(meanx,1/(100*psd)^2)
	      delta~dnorm(0,1/(100*2^0.5*psd)^2)#assuming a common variance the variance of the difference is 2x common var -> 2^0.5 common sd
	      vMinus~dexp(1/29)
	      v<-vMinus+1#v>=1
	}'
	}#}}}
	return (result)
}#}}}


#1= Lee+Wagenmakers using effect size p125
#2= mux versus mux+delta, sigma dunif based on pooled sd, normal(0 or meanx),100*psd or 100*2^0.5*psd) on mux, delta respectively
#3= half cauchy(0,1) on sigma, cauchy(0 or meanx,1) on mux, delta respectively
#4= mux versus mux+delta, sigma dunif based on pooled sd, normal(0 or meanx),100*psd or 100*2^0.5*psd) on mux, delta respectively, robust t-dist on data, Kruschke p462
model<-4
modelString<-getModel(model)
writeLines(modelString,con=modelFile)	
myinits<-switch(model,model1Inits,model2Inits,model3Inits,model4Inits)
data<-switch(model,model1Data,model2Data,model3Data,model2Data)#NB model 4 has same data as model 2
parameters<-switch(model,model1Parameters,model2Parameters,model2Parameters,model4Parameters)#NB models 3 has same parameters as model 2

#model 4 robust
samples4 <- jags(data, inits=myinits, parameters,
	 			model.file =modelFile,
	 			n.chains=3, n.iter=50000, n.burnin=5000, n.thin=1, DIC=T)

model<-2
modelString<-getModel(model)
writeLines(modelString,con=modelFile)	
myinits<-switch(model,model1Inits,model2Inits,model3Inits,model4Inits)
data<-switch(model,model1Data,model2Data,model3Data,model2Data)#NB model 4 has same data as model 2
parameters<-switch(model,model1Parameters,model2Parameters,model2Parameters,model4Parameters)#NB models 3 has same parameters as model 2

#model 2 standard
samples2 <- jags(data, inits=myinits, parameters,
	 			model.file =modelFile,
	 			n.chains=3, n.iter=50000, n.burnin=5000, n.thin=1, DIC=T)
model<-1
modelString<-getModel(model)
writeLines(modelString,con=modelFile)	
myinits<-switch(model,model1Inits,model2Inits,model3Inits,model4Inits)
data<-switch(model,model1Data,model2Data,model3Data,model2Data)#NB model 4 has same data as model 2
parameters<-switch(model,model1Parameters,model2Parameters,model2Parameters,model4Parameters)#NB models 3 has same parameters as model 2

#model 1 effect size
samples1 <- jags(data, inits=myinits, parameters,
	 			model.file =modelFile,
	 			n.chains=3, n.iter=50000, n.burnin=5000, n.thin=1, DIC=T)

save(list = ls(all.names = TRUE), file = session, envir = .GlobalEnv)
}#create session


## ----resultsBayes
library(polspline)
library(BayesFactor)
library(R2jags)
library(BEST)

samples<-ttestBF(dx,dy,posterior=FALSE, iterations=10000)
bftt<-as.data.frame(samples)$bf

#makePlot sets the values of the following, check back after function call
lo95<-0
hi95<-0
noPost<-0
deltaPost<-0
nuPost<-0
rl<-0
rh<-0
sigma<-0
bf<-0

    ciLabel<-function(x0,x1,y,ty,txt){
	lines(c(x0,x1),c(y,y))
	lines(c(x0,x0),c(y-ty,y+ty))
	lines(c(x1,x1),c(y-ty,y+ty))
	text(x0+(x1-x0)/2,y+2*ty, txt)
    }
    multLabel<-function(y0,y1,x,tx,mtx=1){
	arrows(x,y0,x,y1,code=3,length=0.125)
	lines(c(x-2*tx,x),c(y0,y0))
	lines(c(x-2*tx,x),c(y1,y1))
	bf<<-y1/y0
	text(x-mtx*tx,y0+(y1-y0)/2,paste0('x',fn$fndp(bf,2)))
    }


makePlot<-function(samples){#{{{
    # Collect posterior samples across all chains:
    delta.posterior <- samples$BUGSoutput$sims.list$delta  
    mux.posterior <- samples$BUGSoutput$sims.list$mux
    sigma.posterior <- samples$BUGSoutput$sims.list$sigma
    sigma<<-mean(sigma.posterior)
    noPost<<-mean(mux.posterior)
    deltaPost<<-mean(delta.posterior)
    nu.posterior <- NULL
    if('vMinus'%in%samples$parameters.to.save==TRUE){
	nu.posterior <- samples$BUGSoutput$sims.list$vMinus+1
	nuPost<<-mean(nu.posterior)
    }

    fit.posteriorI <- logspline(mux.posterior+mean(delta.posterior))
    fit.posteriorN <- logspline(mux.posterior)

    # 95% confidence interval:
    x0 <- qlogspline(0.025,fit.posteriorI)
    lo95<<-x0
    x1 <- qlogspline(0.975,fit.posteriorI)
    hi95<<-x1
    r0 <- mean(mux.posterior-0.25)#ropeH psd is 2.5 therefore a ROPE covering 0.1 SD around the mean corresponds to small effect
    rl<<-r0
    r1 <- mean(mux.posterior+0.25)#ropeL
    rh<<-r1

    #============ Plot Posteriors ===========================
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
	font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
    xlow  <- -1
    xhigh <- 10.5
    yhigh <- 1.25
    plot(fit.posteriorI, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F, xlab="Rating", ylab="Density")
    axis(1, at = c(-1,0,2,4,6,8,10), lab=c('','0', '2', '4', '6','8', '10'))
    axis(2)
    #mtext(expression(delta), side=1, line = 2.8, cex=2)
    points(jitter(dy), rep(0.025,length(dy)),pch=4, cex=1)
    par(new=T)
    plot (fit.posteriorN, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=2, ylab=" ", xlab = " ") 
    points(jitter(dx), rep(-0.025,length(dx)),pch=3, cex=1)
    legend(x='topright',c('Inh','NoInh'),inset=0.075,lty=1:2, lwd=1:2,pch=c(4,3))

    yv<-dlogspline(mean(mux.posterior),fit.posteriorN)#y value for labels
    ciLabel(x0,x1,yv+0.1,0.035, ' 95% CI')
    ciLabel(r0,r1,yv+0.1,0.035, 'ROPE,')
    return(recordPlot())
}#}}}
pltBayes2<-makePlot(samples2)
lo95_2<-lo95
hi95_2<-hi95
noPost_2<-noPost
deltaPost_2<-deltaPost
sigma_2<-sigma
rl_2<-rl
rh_2<-rh

pltBayes4<-makePlot(samples4)
lo95_4<-lo95
hi95_4<-hi95
noPost_4<-noPost
deltaPost_4<-deltaPost
nuPost_4<-nuPost
sigma_4<-sigma
rl_4<-rl
rh_4<-rh

makePlotBF0<-function(samples){#{{{ plot for Bayes factor
    # Collect posterior samples across all chains:
    delta.posterior <- samples$BUGSoutput$sims.list$delta
    sigma.posterior <- samples$BUGSoutput$sims.list$sigma
    fit.posteriorD <- logspline(delta.posterior)
    deltaPost<<-mean(delta.posterior)
    sigma<<-mean(sigma.posterior)
    # 95% confidence interval:
    x0 <- qlogspline(0.025,fit.posteriorD)
    lo95<<-x0
    x1 <- qlogspline(0.975,fit.posteriorD)
    hi95<<-x1

    #============ Plot Posteriors ===========================
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
	font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
    xlow  <- -1.5
    xhigh <- 2
    yhigh <- 2.5
    plot(fit.posteriorD, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F, xlab=expression(delta), ylab="Density")
    axis(1, at = c(-1.5,-1,0,1,2), lab=c('','-1','0', '1', '2'))
    axis(2)
    par(new=T)
    plot (function(x){dcauchy(x,0,1)}, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=2, ylab=" ", xlab = " ",axes=F) 
    legend(x='topright',c('Posterior','Prior'),inset=0.075,lty=1:2, lwd=1:2)

    yv<-dlogspline(mean(delta.posterior),fit.posteriorD)#y value for labels
    ciLabel(x0,x1,yv+0.1,0.035, ' 95% CI')
    multLabel(dlogspline(0,fit.posteriorD),dcauchy(0),0,0.25,1.4)
    return(recordPlot())
}#}}}
pltBayes1<-makePlotBF0(samples1)
lo95_1<-lo95
hi95_1<-hi95
deltaPost_1<-deltaPost
sigma_1<-sigma

## ----pltB2
pltBayes2

## ----pltB4
pltBayes4

## ----pltB1
pltBayes1


## ----resultsBayes
library(polspline)
library(BayesFactor)
library(R2jags)
library(BEST)

samples<-ttestBF(dx,dy,posterior=FALSE, iterations=10000)
bftt<-as.data.frame(samples)$bf

#makePlot sets the values of the following, check back after function call
lo95<-0
hi95<-0
noPost<-0
deltaPost<-0
nuPost<-0
rl<-0
rh<-0
sigma<-0
bf<-0

    ciLabel<-function(x0,x1,y,ty,txt){
	lines(c(x0,x1),c(y,y))
	lines(c(x0,x0),c(y-ty,y+ty))
	lines(c(x1,x1),c(y-ty,y+ty))
	text(x0+(x1-x0)/2,y+2*ty, txt)
    }
    multLabel<-function(y0,y1,x,tx,mtx=1){
	arrows(x,y0,x,y1,code=3,length=0.125)
	lines(c(x-2*tx,x),c(y0,y0))
	lines(c(x-2*tx,x),c(y1,y1))
	bf<<-y1/y0
	text(x-mtx*tx,y0+(y1-y0)/2,paste0('x',fn$fndp(bf,2)))
    }


makePlot<-function(samples){#{{{
    # Collect posterior samples across all chains:
    delta.posterior <- samples$BUGSoutput$sims.list$delta  
    mux.posterior <- samples$BUGSoutput$sims.list$mux
    sigma.posterior <- samples$BUGSoutput$sims.list$sigma
    sigma<<-mean(sigma.posterior)
    noPost<<-mean(mux.posterior)
    deltaPost<<-mean(delta.posterior)
    nu.posterior <- NULL
    if('vMinus'%in%samples$parameters.to.save==TRUE){
	nu.posterior <- samples$BUGSoutput$sims.list$vMinus+1
	nuPost<<-mean(nu.posterior)
    }

    fit.posteriorI <- logspline(mux.posterior+mean(delta.posterior))
    fit.posteriorN <- logspline(mux.posterior)

    # 95% confidence interval:
    x0 <- qlogspline(0.025,fit.posteriorI)
    lo95<<-x0
    x1 <- qlogspline(0.975,fit.posteriorI)
    hi95<<-x1
    r0 <- mean(mux.posterior-0.25)#ropeH psd is 2.5 therefore a ROPE covering 0.1 SD around the mean corresponds to small effect
    rl<<-r0
    r1 <- mean(mux.posterior+0.25)#ropeL
    rh<<-r1

    #============ Plot Posteriors ===========================
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
	font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
    xlow  <- -1
    xhigh <- 10.5
    yhigh <- 1.25
    plot(fit.posteriorI, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F, xlab="Rating", ylab="Density")
    axis(1, at = c(-1,0,2,4,6,8,10), lab=c('','0', '2', '4', '6','8', '10'))
    axis(2)
    #mtext(expression(delta), side=1, line = 2.8, cex=2)
    points(jitter(dy), rep(0.025,length(dy)),pch=4, cex=1)
    par(new=T)
    plot (fit.posteriorN, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=2, ylab=" ", xlab = " ") 
    points(jitter(dx), rep(-0.025,length(dx)),pch=3, cex=1)
    legend(x='topright',c('Inh','NoInh'),inset=0.075,lty=1:2, lwd=1:2,pch=c(4,3))

    yv<-dlogspline(mean(mux.posterior),fit.posteriorN)#y value for labels
    ciLabel(x0,x1,yv+0.1,0.035, ' 95% CI')
    ciLabel(r0,r1,yv+0.1,0.035, 'ROPE,')
    return(recordPlot())
}#}}}
pltBayes2<-makePlot(samples2)
lo95_2<-lo95
hi95_2<-hi95
noPost_2<-noPost
deltaPost_2<-deltaPost
sigma_2<-sigma
rl_2<-rl
rh_2<-rh

pltBayes4<-makePlot(samples4)
lo95_4<-lo95
hi95_4<-hi95
noPost_4<-noPost
deltaPost_4<-deltaPost
nuPost_4<-nuPost
sigma_4<-sigma
rl_4<-rl
rh_4<-rh

makePlotBF0<-function(samples){#{{{ plot for Bayes factor
    # Collect posterior samples across all chains:
    delta.posterior <- samples$BUGSoutput$sims.list$delta
    sigma.posterior <- samples$BUGSoutput$sims.list$sigma
    fit.posteriorD <- logspline(delta.posterior)
    deltaPost<<-mean(delta.posterior)
    sigma<<-mean(sigma.posterior)
    # 95% confidence interval:
    x0 <- qlogspline(0.025,fit.posteriorD)
    lo95<<-x0
    x1 <- qlogspline(0.975,fit.posteriorD)
    hi95<<-x1

    #============ Plot Posteriors ===========================
    par(cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5,
	font.lab = 2, cex.axis = 1.3, bty = "n", las=1)
    xlow  <- -1.5
    xhigh <- 2
    yhigh <- 2.5
    plot(fit.posteriorD, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lty=1, lwd=1, axes=F, xlab=expression(delta), ylab="Density")
    axis(1, at = c(-1.5,-1,0,1,2), lab=c('','-1','0', '1', '2'))
    axis(2)
    par(new=T)
    plot (function(x){dcauchy(x,0,1)}, ylim=c(0,yhigh), xlim=c(xlow,xhigh), lwd=2, lty=2, ylab=" ", xlab = " ",axes=F) 
    legend(x='topright',c('Posterior','Prior'),inset=0.075,lty=1:2, lwd=1:2)

    yv<-dlogspline(mean(delta.posterior),fit.posteriorD)#y value for labels
    ciLabel(x0,x1,yv+0.1,0.035, ' 95% CI')
    multLabel(dlogspline(0,fit.posteriorD),dcauchy(0),0,0.25,1.4)
    return(recordPlot())
}#}}}
pltBayes1<-makePlotBF0(samples1)
lo95_1<-lo95
hi95_1<-hi95
deltaPost_1<-deltaPost
sigma_1<-sigma

## ----pltB2
pltBayes2

## ----pltB4
pltBayes4

## ----pltB1
pltBayes1

## ----extDataHDI
#look at posterior distn for the difference in responding during extinction for inhibitors and non-inhibitors
ext<-new.env()
evalq({
restore<-TRUE
load("extData.RData")
resps<-with(extData,aggregate(xrsp1,list(sno,igrp),sum))
colnames(resps)<-c('sno','igrp','nx')
nxi<-resps$nx[resps$igrp==1]
meannxI<-mean(nxi)
nxn<-resps$nx[resps$igrp==0]
nex<-8#n extinction trials
ni<-length(nxi)#n inhibitors
nn<-length(nxn)#n non-inhibitors
#pooled sd 
psd<-(((ni-1)*var(nxi)+(nn-1)*var(nxn))/(ni+nn-2))^0.5

session<-'extSamples.RData'
if(restore==TRUE)
    load(session)
else{#only do the MCMC if needed
    modelStr<-'model{
	for(i in 1:ni){
		nxi[i]~dbin(thetaI[i],nex)
		thetaI[i]<-phi(phiI[i])
		phiI[i] ~ dnorm(muI,lambda)
	    }
	for(j in 1:nn){
		nxn[j]~dbin(thetaN[j],nex)
		thetaN[j]<-phi(phiN[j])
		phiN[j] ~ dnorm(muI+delta,lambda)
	    }
	#priors
	sigma~dunif(psd/100,psd*100)
	lambda<-pow(sigma,-2)
	muI~dnorm(meannxI,1/(100*psd)^2)
	delta~dnorm(0,1/(100*2^0.5*psd)^2)#assuming a common variance the variance of the difference is 2x common var -> 2^0.5 common sd
    }'
    writeLines(modelStr,con=modelFile)	
    data<-list('nxi','nxn','nex','ni','nn','psd','meannxI')
    inits<-list(
      list(delta = rnorm(1,0,2), muI = rnorm(1,0.2,2), sigma = runif(1,0,5)),
      list(delta = rnorm(1,0,2), muI = rnorm(1,0.2,2), sigma = runif(1,0,5)),
      list(delta = rnorm(1,0,2), muI = rnorm(1,0.2,2), sigma = runif(1,0,5)))
    params<-c('delta','muI','sigma')

    extSamples<-jags(data,inits=inits,params,model.file=modelFile,n.chains=3,n.iter=10000,n.burnin=1000,n.thin=1,DIC=T)
    save(extSamples,file=session)
}#create session
delta.posterior <- extSamples$BUGSoutput$sims.list$delta      
muI.posterior <- extSamples$BUGSoutput$sims.list$muI      
sig.posterior <- extSamples$BUGSoutput$sims.list$sigma
    
#plot(logspline(muI.posterior),xlim=c(-1.25,0))
#plot(logspline(muI.posterior+mean(delta.posterior)),add=TRUE)

fit.posteriorN <- logspline(muI.posterior+mean(delta.posterior))

    # 95% confidence interval:
    x0 <- qlogspline(0.025,fit.posteriorN)
    x0_p<-pnorm(x0,0,1)#map to prop
    x1 <- qlogspline(0.975,fit.posteriorN)
    x1_p <- pnorm(x1,0,1)#map to prop
    r0 <- mean(muI.posterior)-0.1*mean(sig.posterior)#ropeL a ROPE covering 0.1 SD around the mean corresponds to small effect
    r0_p<-pnorm(r0,0,1)
    r1<- mean(muI.posterior)+0.1*mean(sig.posterior)#ropeH
    r1_p<-pnorm(r1,0,1)
},envir=ext)
