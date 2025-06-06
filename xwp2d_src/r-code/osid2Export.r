## ----setup
rm(list=ls(all=TRUE))
setDirs<-function(){#{{{
    cwd<-getwd()
    rstr<-'c:\\users\\admin\\'
    #rstr<-'c:\\users\\spg\\'
    setwd(paste0(rstr,'documents\\rwork\\e168'))
    source('..\\functions\\generic.r')#assumes we are currently in rwork\XXX
    library(Hmisc) 
    return(cwd)
}#}}}
cwd<-setDirs()
#functions and variables accessed by e1 and e2
symbols<-c(15,0,16,1,17,2,18,5) 
ymin<-0#was ymin<--0.5 why?
ymax<-1 
yMax<-10.5 
gcex<-1.1 

#Start data functions#{{{
getP1<-function(df,c){
  result<-with(df[df$cue==c,], aggregate(xrsp1,list(sno,block),mean))
  colnames(result)<-c('sno','block','xrsp')
  return(result)
}

getMeansSE<-function(df){
    result<-with(df, aggregate(xrsp, list(block),mean))
    colnames(result)<-c('block','xrsp')
    tmp<-with(df, aggregate(xrsp, list(block),serr))
    result$se<-tmp$x
    return(result)
}

#return inhibitor IDs
#inhibitor for those who have nr<n responses in stage 2a G tests
#nr[0,1,2] so n should be 1 or 2
inhn<-1
getInh<-function(df,n){
    df<-with(df[df$block==6|df$block==7,],aggregate(xrsp,list(sno),sum))
    colnames(df)<-c('sno','nr')
    df$inh<-0
    df$inh[df$nr<n]<-1
    return(df$sno[df$inh==1])
}


getFNData<-function(df){
    S1<-with(df, aggregate(df$rsp,list(sno,block,cue),mean))
    colnames(S1)<-c('sno','block','cue','rsp')
    means<-with(S1, aggregate(rsp, list(block,cue),mean))
    colnames(means)<-c('block','cue','rsp')
    df<-with(S1, aggregate(rsp, list(block,cue),serr))
    means$se<-df$x
    return(means)
}

getInhData<-function(df){
    df<-df[df$block==11|df$block==12,]
    df<-with(df,aggregate(rsp,list(sno),mean))
    colnames(df)<-c('sno','rsp')
    return(df)
}
#End data functions#}}}

#Start plotting functions#{{{
plotTrace<-function(df, xl, yl){
    plt<-xyplot(xrsp~block,data=df, type='b', 
		ylab=list(label=yl,cex=gcex), 
	        xlab=xl,ylim=c(ymin,ymax), se=df$se, panel=fn$mypanel)
}

plotFN<-function(df,xl,yl){
    means<-getFNData(df)
    return(xyplot(rsp~block, groups=cue,data=means, type='b', 
		    ylab=list(label=xl,cex=gcex), 
		    xlab=list(label=yl,cex=gcex),
		    ylim=c(ymin,yMax), se=means$se, 
		    panel=fn$mypanelG,
		    key=list(space='right',type=c('p'), 
				  text=list(c('A+','AB-','C-','CD+')),
				  points=list(pch=symbols[1:4]))))
}

plotFNReverse<-function(df,xl,yl){
    means<-getFNData(df)
    return(xyplot(rsp~block, groups=cue,data=means, type='b', 
		    ylab=list(label=xl,cex=gcex), 
		    xlab=list(label=yl,cex=gcex),
		    ylim=c(ymin,yMax), se=means$se, 
		    panel=fn$mypanelG,
		    key=list(space='right',
				  text=list(c('B+','M-')),
				  points=list(pch=symbols[1:2])),
		                  scales = list(x = list(at = c(5:8),labels=c(5:8)))))
}

plotFNReverseGrp<-function(df,xl,yl,ids){
    means1<-getFNData(df[df$sno%in%ids,])
    means1$igrp<-1
    means0<-getFNData(df[!df$sno%in%ids,])
    means0$igrp<-0
    means<-rbind(means0,means1)
    means$igrp<-factor(means$igrp,labels=c('NoInh','Inh'))
    return(xyplot(rsp~block|igrp, groups=cue,data=means, type='b', 
		    ylab=list(label=xl,cex=gcex), 
		    xlab=list(label=yl,cex=gcex),
		    ylim=c(ymin,yMax), se=means$se, 
		    panel=fn$mypanelG,
		    key=list(space='right',
	  	    text=list(c('B+','M-')),
		    points=list(pch=symbols[1:2])),
		    scales = list(x = list(at = c(5:8),labels=c(5:8)))))
}

plotFNTest<-function(df,xl,yl){
    means<-getFNData(df)
    return(xyplot(rsp~block,data=means, type='b',
	      ylab=list(label=yl,cex=gcex), 
	      xlab=list(label=xl,cex=gcex),
	      ylim=c(ymin,yMax), se=means$se, panel=fn$mypanel,
	      scales = list(x = list(at = c(9:12),labels=c('A+','A+','AB-','AB-')))))
}

plotInhTest<-function(df,xl,yl){
    means<-with(df, aggregate(rsp, list(igrp),mean))
    colnames(means)<-c('grp','rsp')
    tmp<-with(df, aggregate(rsp, list(igrp),serr))
    means$se<-tmp$x
    means$grp<-as.factor(means$grp)
    return(barchart(rsp~grp, ylab=yl, xlab=xl, se=means$se, xlim=c(0.5,2.5),ylim=c(0,3.5), 
		    panel=fn$mybar,data=means,scales = list(x = list(at = c(1,2),labels=c('NoInh','Inh')))))
}

plotNP<-function(df,xl,yl){
    means<-getFNData(df)
    return(xyplot(rsp~block, groups=cue,data=means, type='b', 
		    ylab=list(label=xl,cex=gcex), 
		    xlab=list(label=yl,cex=gcex),
		    ylim=c(ymin,yMax), se=means$se, 
		    panel=fn$mypanelG,
		    scales = list(x = list(at = c(13:15),labels=c(13:15))),
		    key=list(space='right',type=c('p'), 
				  text=list(c('A+','AB-','B+')),
				  points=list(pch=symbols[1:3]))))
}

plotNPGrp<-function(df,xl,yl,ids){
    means1<-getFNData(df[df$sno%in%ids,])
    means1$igrp<-1
    means0<-getFNData(df[!df$sno%in%ids,])
    means0$igrp<-0
    means<-rbind(means0,means1)
    means$igrp<-factor(means$igrp,labels=c('NoInh','Inh'))
    return(xyplot(rsp~block|igrp, groups=cue,data=means, type='b', 
		    ylab=list(label=xl,cex=gcex), 
		    xlab=list(label=yl,cex=gcex),
		    ylim=c(ymin,yMax), se=means$se, 
		    panel=fn$mypanelG,
		    scales = list(x = list(at = c(13:15),labels=c(13:15))),
		    key=list(space='right',type=c('p'), 
				  text=list(c('A+','AB-','B+')),
				  points=list(pch=symbols[1:3]))))
}

#plotHists helper
ph<-function(df){
    br<-c(0,2,4,6,8,10)
    hist(df$rsp,breaks=br)
    hist(df$rsp[df$igrp==0],breaks=br)
    hist(df$rsp[df$igrp==1],breaks=br)
}
plotHists<-function(){
    par(mfcol=c(3,2))
    ph(e1$idata)
    ph(e2$idata)
    result<-recordPlot()
    plot.new()
    return(result)
}
#End plotting functions#}}}

#Start experiment 1#{{{
e1<-new.env()
evalq({
load('AllData168a.RData') #part 1
load('AllData168b.RData') #part 2
colnames(AllDataA)[3]<-'cue'#make df name compatible with E175
AllDataB$block[AllDataB$block>10]<-AllDataB$block[AllDataB$block>10]+2#make block numbering compatible with E175
AllDataB$block[AllDataB$trial==57]<-10
AllDataB$block[AllDataB$trial==58]<-11
AllDataB$block[AllDataB$trial==59]<-12

#PART 1 
ATrials1<-subset(AllDataA,cue=='A') 
#recode the recovery test block value 
ATrials1$block[ATrials1$block==11 & ATrials1$twb==0]<-10 
AggA1<-getP1(ATrials1,'A')
MeansA1<-getMeansSE(AggA1)
TAPlot<-plotTrace(MeansA1,'','Proportion X responses')

GTrials1<-subset(AllDataA,cue=='G') 
#recode test trial block 
GTrials1$block[GTrials1$block==10]<-7 
GTrials1$block[GTrials1$block==7 & GTrials1$twb==0]<-6
AggG1<-getP1(GTrials1,'G')
iids<-getInh(AggG1,inhn)
MeansG1<-getMeansSE(AggG1)
TGPlot<-plotTrace(MeansG1,'','Proportion X responses')

#PART 2
FNPlot<-plotFN(AllDataB[AllDataB$block<5,],'Rating','')
FNRPlot<-plotFNReverse(AllDataB[AllDataB$block>4&AllDataB$block<9,],'Rating','')
FNRPlotGrp<-plotFNReverseGrp(AllDataB[AllDataB$block>4&AllDataB$block<9,],'Rating','',iids)
FNTPlot<-plotFNTest(AllDataB[AllDataB$block>8&AllDataB$block<13,],'Rating','')

#inh test
idata<-getInhData(AllDataB)
idata$igrp<-0
idata$igrp[idata$sno%in%iids]<-1
ITestPlot<-plotInhTest(idata,'Group','Rating')
},envir=e1)
#End experiment 1#}}}

setwd('..\\e175')

#Start experiment 2#{{{
e2<-new.env()
evalq({
load('part1.RData')
load('part2.RData')
#load('part3.RData')
colnames(p1)[10]<-'xrsp1'#make df name compatible with E168

#PART 1 
ATrials1<-subset(p1,cue=='A')
#recode the recovery test block value 
ATrials1$block[ATrials1$block==12]<-10
ATrials1$block[ATrials1$block==13]<-11
AggA1<-getP1(ATrials1,'A')
MeansA1<-getMeansSE(AggA1)
TAPlot<-plotTrace(MeansA1,'','Proportion X responses')

GTrials1<-subset(p1,cue=='G') 
#recode test trial block 
GTrials1$block[GTrials1$block==10]<-6
GTrials1$block[GTrials1$block==11]<-7
AggG1<-getP1(GTrials1,'G')
iids<-getInh(AggG1,inhn)
MeansG1<-getMeansSE(AggG1)
TGPlot<-plotTrace(MeansG1,'','Proportion X responses')

#PART 2
FNPlot<-plotFN(p2[p2$block<5,],'Rating','')
FNRPlot<-plotFNReverse(p2[p2$block>4&p2$block<9,],'Rating','')
FNRPlotGrp<-plotFNReverseGrp(p2[p2$block>4&p2$block<9,],'Rating','',iids)
FNTPlot<-plotFNTest(p2[p2$block>8&p2$block<13,],'Rating','')

#inh test
idata<-getInhData(p2)
idata$igrp<-0
idata$igrp[idata$sno%in%iids]<-1
ITestPlot<-plotInhTest(idata,'Group','Rating')
},envir=e2)
#End experiment 2#}}}

writeData<-FALSE
writeITestData<-function(fn){
    ex1<-e1$idata
    ex2<-e2$idata
    ex2$sno<-ex2$sno+100
    write.table(rbind(ex1,ex2),fn,row.names=F,col.names=F,sep='\t')	
}
if(writeData==TRUE)
    writeITestData('itest.txt')

#glautierEtAl2013
#getting the means for the combined extinction groups and combined controls for comparison with current investigation
library(foreign)
#.sav data from adminpc F:\d\ppexpts\expt133\all.sav, Experiment 1, see also reade133.sps from that directory
E1prev<-read.spss(file='../../latex/work/osid/allE133.sav',to.data.frame=TRUE)
#with(E1prev,aggregate(avG6,list(grp),mean))
#with(E1prev,aggregate(avG6,list(grp),serr))
#with(E1prev,aggregate(avG6,list(grp),length))
#see E147/E147#1.r and E147#2.r
load('../E147/AllData.RData')
#with(AggG[AggG$block==10,],aggregate(pxrsp,list(grp),mean))
#with(AggG[AggG$block==10,],aggregate(pxrsp,list(grp),serr))
#with(AggG[AggG$block==10,],aggregate(pxrsp,list(grp),length))
control<-c(E1prev$avG6[E1prev$grp=='noext'],with(AggG[AggG$block==10,],pxrsp[grp=='g0']),with(AggG[AggG$block==10,],pxrsp[grp=='g3']))
ext<-c(E1prev$avG6[E1prev$grp=='sext'],E1prev$avG6[E1prev$grp=='mext'],
       with(AggG[AggG$block==10,],pxrsp[grp=='g1']),with(AggG[AggG$block==10,],pxrsp[grp=='g2']))
gv<-function(v){
    c(mean(v),serr(v),length(v))
}

ext1<-c(with(e1$AggG1[e1$AggG1$block>5,],aggregate(xrsp,list(sno),mean))$x,
	with(e2$AggG1[e2$AggG1$block>5,],aggregate(xrsp,list(sno),mean))$x)
ciDataPrev<-gv(control)
eiDataPrev<-gv(ext)
eiDataCurr<-gv(ext1)
setwd(cwd)

## ----TRACEAoE1
e1$TAPlot

## ----TRACEGoE1
e1$TGPlot

## ----FNE1
e1$FNPlot

## ----FNRE1
e1$FNRPlot

## ----FNRGE1
e1$FNRPlotGrp

## ----FNTE1
e1$FNTPlot

## ----InTestE1
e1$ITestPlot

## ----TRACEAoE2
e2$TAPlot

## ----TRACEGoE2
e2$TGPlot

## ----FNE2
e2$FNPlot

## ----FNRE2
e2$FNRPlot

## ----FNRGE2
e2$FNRPlotGrp

## ----FNTE2
e2$FNTPlot

## ----InTestE2
e2$ITestPlot

## ----InTestsAll
idataAll<-rbind(e1$idata,e2$idata)
plotInhTest(idataAll,'Group','Rating')

## ----inddiff{{{
#inhibitors also show stronger extinction
getID<-function(d1,d2){
    t1<-d1
    t1$inh<-0
    t1$inh[t1$sno%in%e1$iids]<-1
    t2<-d2
    t2$inh<-0
    t2$inh[t2$sno%in%e2$iids]<-1
    t2$sno<-t2$sno+100
    ixr<-rbind(t1[,c('sno','block','xrsp1','inh')],t2[,c('sno','block','xrsp1','inh')])
    pind<-with(ixr,aggregate(xrsp1,list(block,inh),mean))
    pind$se<-with(ixr,aggregate(xrsp1,list(block,inh),serr))$x
    colnames(pind)<-c('block','grp','xrsp','se')
    return(pind)
}

plotID<-function(d,t){
    xyplot(xrsp~block,groups=grp,data=d,main=t,type='b',pch=symbols,key=list(space='right',type=c('p'),text=list(c('NoInh','Inh')),points=list(pch=symbols[1:2]))) 
}

AData<-getID(e1$ATrials1,e2$ATrials1)
GData<-getID(e1$GTrials1,e2$GTrials1)
plotID(AData,'ATrials')
plotID(GData,'GTrials')
}}}
