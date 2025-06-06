#create an environment and store in it a bunch of useful functions
library(lattice)
library(reshape2)
library(ez)
library(xtable)
library('expm')

if('fn' %in% search())
	detach(fn)
fn<-new.env()

fn$first<-function(x){return (head(x,1))}
fn$serr<-function(x){sqrt(var(x)/length(x))}

fn$mypanel<-function(x, y, se, subscripts, pch=15, cex=1.5,...){
 x<-as.numeric(x)
 y<-as.numeric(y)
 ly<-y-se[subscripts]
 uy<-y+se[subscripts]
 panel.arrows(x,ly,x,uy, col='black', length=0.1, angle=90, code=3)
 panel.xyplot(x,y,type='b',pch=pch, col='black',cex=cex)
 #print(se[subscripts])
}

fn$mypanelG<-function(x, y, se, groups, subscripts, ...){
	x<-as.numeric(x)
	y<-as.numeric(y)
	ly<-y-se[subscripts]
	uy<-y+se[subscripts]
	panel.arrows(x,ly,x,uy, col='black', length=0.1, angle=90, code=3)
	panel.xyplot(x,y,type='b',pch=c(15,0,16,1,17,2,18,5), col='black', cex=1.5,groups=groups, subscripts=subscripts)
	#print(se[subscripts])
}

fn$mypanelG1<-function(x, y, se, groups, subscripts, ...){
	x<-as.numeric(x)
	y<-as.numeric(y)
	ly<-y-se[subscripts]
	uy<-y+se[subscripts]
	panel.arrows(x,ly,x,uy, col='black', length=0.1, angle=90, code=3)
	panel.xyplot(x,y,type='b',pch=c(15,1,2,0,16,18,5,17), col='black', cex=1.5,groups=groups, subscripts=subscripts)
	#print(se[subscripts])
}

fn$mypanelG2<-function(x, y, se, groups, subscripts, ...){
	x<-as.numeric(x)
	y<-as.numeric(y)
	ly<-y-se[subscripts]
	uy<-y+se[subscripts]
	panel.arrows(x,ly,x,uy, col='black', length=0.1, angle=90, code=3)
	panel.xyplot(x,y,type='b',pch=c(1,2,0,16,18,5,17,15), col='black', cex=1.5,groups=groups, subscripts=subscripts)
	#print(se[subscripts])
}

fn$mybar<-function(x,y,se,subscripts,...){
 uy<-y+se[subscripts]
 panel.arrows(x,y,x,uy, col='black', length=0.1, angle=90, code=3)
 panel.barchart(x,y, horizontal=FALSE, ylim=c(0,6))
} 

#pass eos - error bar offsets of length x (i.e. four bars, four values), adjust eos to fit the plot
fn$mybarG<-function(x,y,se,groups,subscripts,eos,...){
 uy<-y+se[subscripts]
 xl<-as.numeric(x)
 xl<-eos+xl
 panel.arrows(xl,y,xl,uy, col='black', length=0.1, angle=90, code=3)
 panel.barchart(x,y,groups=groups, subscripts=subscripts,horizontal=FALSE)
} 
#format results of t-test
fn$tout<-function(tt){
	#df will have up to one dp and drop trailing zeroes
	#t value will have 3 significant digits
	#p value will have one of four forms
	if(tt$p.value<.05)
		eq<-', p<.05'
	if(tt$p.value<.01)
		eq<-', p<.01'
	if(tt$p.value<.001)
		eq<-', p<.001'
	if(tt$p.value>=.05)
		eq<-paste(', p=',formatC(tt$p.value,format='f', digits=2,drop0trailing=TRUE),sep='')
	result<-paste('$t(',formatC(tt$parameter,format='f',digits=1,drop0trailing=TRUE),')=',format(abs(tt$statistic),digits=3),eq,'$',sep='')
	result
}

#format F value
#e.g. fout(1.784,0.5014483*24,0.5014483*840)
fn$fout<-function(F, dfn, dfd){
	#df will have up to one dp and drop trailing zeroes
	#F value will have 3 significant digits
	#p value will have one of four forms
	p.value<-(1-pf(F,dfn, dfd))
	if(p.value<.05)
		eq<-', p<.05'
	if(p.value<.01)
		eq<-', p<.01'
	if(p.value<.001)
		eq<-', p<.001'
	if(p.value>=.05)
		eq<-paste(', p=',formatC(p.value,format='f', digits=2,drop0trailing=TRUE),sep='')
	result<-paste('$F(',formatC(dfn,format='f',digits=1,drop0trailing=TRUE),',',formatC(dfd,format='f',digits=1,drop0trailing=TRUE),')=',format(abs(F),digits=3),eq,'$',sep='')
	result
}

#format F value with effect size - generalised eta squared
#e.g. fout(1.784,0.5014483*24,0.5014483*840, 0.35)
fn$foute<-function(F, dfn, dfd,g){
	#df will have up to one dp and drop trailing zeroes
	#F value will have 3 significant digits
	#p value will have one of four forms
	p.value<-(1-pf(F,dfn, dfd))
	if(p.value<.05)
		eq<-', p<.05'
	if(p.value<.01)
		eq<-', p<.01'
	if(p.value<.001)
		eq<-', p<.001'
	if(p.value>=.05)
		eq<-paste(', p=',formatC(p.value,format='f', digits=2,drop0trailing=TRUE),sep='')
	if(g>=.005)
		result<-paste('$F(',formatC(dfn,format='f',digits=1,drop0trailing=TRUE),',',formatC(dfd,format='f',digits=1,drop0trailing=TRUE),')=',format(abs(F),digits=3),eq,',\\hat{\\eta}_{G}^{2}=',sub('0*\\.','\\.',formatC(g,format='f', digits=2,drop0trailing=TRUE)),'$',sep='')
	    else
		result<-paste('$F(',formatC(dfn,format='f',digits=1,drop0trailing=TRUE),',',formatC(dfd,format='f',digits=1,drop0trailing=TRUE),')=',format(abs(F),digits=3),eq,',\\hat{\\eta}_{G}^{2}<.01$',sep='')
	result
}

fn$formatp<-function(p){
	if(p<.05)
		eq<-', p<.05'
	if(p<.01)
		eq<-', p<.01'
	if(p<.001)
		eq<-', p<.001'
	if(p>=.05)
		eq<-paste(', p=',formatC(p,format='f', digits=2,drop0trailing=TRUE),sep='')
	return(eq)
}

#request statistic s, for effect e, from detailed ezANOVA a
#result<-ezANOVA(data=data1,dv=.(pc),wid=.(sno),within=.(g,b,s),between=.(grp),detailed=TRUE,return_aov=TRUE)
#e.g. getAI('SSd','g:b',result) will get SSd for the g x b interaction from result
fn$getAI<-function(s,e,a){
	return(a[[1]][[s]][match(e,a[[1]][['Effect']])])
}

#given a detailed ezANOVA provide formatted f statistic e.g. getFE('igrp:block',result)
fn$getFE<-function(e,a){
    F<-getAI('F',e,a)
    DFn<-getAI('DFn',e,a)
    DFd<-getAI('DFd',e,a)
    g<-getAI('ges',e,a)
    idx<-match(e,a[[2]][['Effect']])
    if(is.na(idx)==TRUE)
       	return(foute(F,DFn,DFd,g))#return without sphericity correction, the effect does not appear in sphericity tests
    pgge<-a[[2]][['p']][idx]
    if(pgge>=0.05)
       	return(foute(F,DFn,DFd,g))#return without sphericity correction, sphericity is not violated
    gge<-a[[3]][['GGe']][idx]
    return(foute(F,DFn*gge,DFd*gge,g))#return with sphericity correction
}
#given a detailed ezANOVA provide formatted f statistic e.g. getFE('igrp:block',result)
fn$getFENoSphericity<-function(e,a){
    F<-getAI('F',e,a)
    DFn<-getAI('DFn',e,a)
    DFd<-getAI('DFd',e,a)
    g<-getAI('ges',e,a)
    return(foute(F,DFn,DFd,g))
}

#for df created by object of class PegsData or PegsDataSel
#NB create response and progO in df before passing, see readE166.r
fn$getFlip<-function(data){
	ids<-unique(data$sno)
	n<-length(ids)
	result<-data[0,]
	for(i in 1:n){
		buff<-data[data$sno==ids[i],]
		buff$flip<-0
	 	buff$flip[buff$progO!=buff$outcomeOrder]<-1
		if(1 %in% buff$flip){#flip Ys and Xs
			buff$response[buff$response=='Y']<-'X1'
			buff$response[buff$response=='X']<-'Y'
			buff$response[buff$response=='X1']<-'X'
		}
		result<-rbind(result,buff)
	}
	result
}

#for df created by object of class PegsData or PegsDataSel
#this only works for experiments where single outcome X is programmed
#fn$getFlip<-function(data){
#	ids<-unique(data$sno)
#	n<-length(ids)
#	result<-data[0,]
#	for(i in 1:n){
#		buff<-data[data$sno==ids[i],]
#		if('Y' %in% buff$outcomeOrder){#flip Ys and Xs
#			buff$response[buff$response=='Y']<-'X1'
#			buff$response[buff$response=='X']<-'Y'
#			buff$response[buff$response=='X1']<-'X'
#		}
#		result<-rbind(result,buff)
#	}
#	result
#}

#a random r x c matrix, all values min<=v<=max drawn from uniform distribution
#r,c integers > 0
fn$randMatrix<-function(r=3,c=3,min=0, max=1){
  matrix(c(runif(r*c)),r,c)
}

#a random r x c matrix all values 0<=v<=1 and all row sums=1
#r,c integers > 0
fn$randMatrixRSum1<-function(r=3,c=3){
if(c==1)
	return(matrix(data=1,nrow=r,ncol=c))
getRow<-function(c){
  row<-vector(mode='numeric', length=c)
  for(i in 1:(c-1)){
    if(i==1)
		row[i]<-runif(1)
	else{
		row[i]<-runif(1,
					max=(function(i){
						result<-0
						for(j in 1:(i-1)){
							result<-result+row[j]
						}
						return(1-result)
					})(i))
		}
	}
	row[c]<-(function(){
		result<-0
		for(i in 1:(c-1))
			result<-result+row[i]
		return(1-result)})()
	return(sample(row))#randomise order
}
result<-matrix(data=NA,nrow=r,ncol=c)
for(i in 1:r)
	result[i,]<-getRow(c)
return(result)
}

#function first returns the first row of the dataframe passed
fn$first<-function(x){return (head(x,1))}

#return the lastN characters of a string
fn$lastN<-function(string, n){
	l<-nchar(string)
	result<-substr(string,l-n+1,l)
	result
}

#return first characters of a string, chopping-off last n
fn$loseLastN<-function(string,n){
	l<-nchar(string)
	result<-substr(string,0,l-n)
	result
}


#http://wiki.cbr.washington.edu/qerm/index.php/R/Contour_Plots
#http://wiki.cbr.washington.edu/qerm/sites/qerm/images/4/44/Filled.contour2.R
#http://stackoverflow.com/questions/16812528/filled-contour-in-r-3-0-x-throws-error
#filled.contour2
fn$filled.contour2 <-#{{{
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2]) * par("csi") * 2.54
  par(las = las)
  mar <- mar.orig
  plot.new()
  par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}#}}}
#http://wiki.cbr.washington.edu/qerm/sites/qerm/images/1/16/Filled.contour3.R
#filled.contour3
fn$filled.contour3 <-#{{{
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
{
  # modification by Ian Taylor of the filled.contour function
  # to remove the key and facilitate overplotting with contour()
  # further modified by Carey McGilliard and Bridget Ferris
  # to allow multiple plots on one page

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
 # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
 # on.exit(par(par.orig))
 # w <- (3 + mar.orig[2]) * par("csi") * 2.54
 # par(las = las)
 # mar <- mar.orig
 plot.new()
 # par(mar=mar)
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
    stop("no proper 'z' matrix specified")
  if (!is.double(z)) 
    storage.mode(z) <- "double"
  .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                          col = col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}#}}}

#http://wiki.cbr.washington.edu/qerm/sites/qerm/images/2/25/Filled.legend.R
#filled.legend
fn$filled.legend <-#{{{
  function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
    length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) 
{
  # modification of filled.contour by Carey McGilliard and Bridget Ferris
  # designed to just plot the legend
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq.int(0, 1, length.out = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
  #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  #  on.exit(par(par.orig))
  #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  #  par(las = las)
  #  mar <- mar.orig
  #  mar[4L] <- mar[2L]
  #  mar[2L] <- 1
  #  par(mar = mar)
   # plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
        yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
}#}}}
    #
#    if (!missing(key.title)) 
#        key.title
#    mar <- mar.orig
#    mar[4L] <- 1
#    par(mar = mar)
#    plot.new()
#    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
#    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
#        stop("no proper 'z' matrix specified")
#    if (!is.double(z)) 
#        storage.mode(z) <- "double"
#    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
#        col = col))
#    if (missing(plot.axes)) {
#        if (axes) {
#            title(main = "", xlab = "", ylab = "")
#            Axis(x, side = 1)
#            Axis(y, side = 2)
#        }
#    }
#    else plot.axes
#    if (frame.plot) 
#        box()
#    if (missing(plot.title)) 
#        title(...)
#    else plot.title
#    invisible()
#}


# patchSynctex(texfile) - for texfile="blubb.tex" it does use the "blubb-concordance.tex" 
#  file generated by knitr to patch the blubb.synctex[.gz] file generated by latex in order to
#  fix forward and reverse seach between blubb.pdf and blubb.Rnw
#
# tested only under windows, with one single utf8 Rnw file.
#
# Author: Jan Gleixner <jan.gleixner@gmail>
###############################################################################

fn$patchKnitrSynctex <- function (texfile){
	require(tools)
	f=paste0(tools::file_path_sans_ext(texfile), "-concordance.tex")
	if (!file.exists(f)) 
		stop(paste(f,"does not exist! Did you set 'opts_knit$set(concordance = TRUE);'?"))
	text<-readChar(f, file.info(f)$size);
	require(stringr)
	re="\\\\Sconcordance\\{concordance:([^:]*):([^\\%]*):\\%\\r?\\n(\\d+)(( \\d+ \\d+)*)\\}";
	parsed=str_match_all(text,re);
	for(i in seq(1,nrow(parsed[[1]]))){		
		texF=parsed[[1]][i,2];
		rnwF=parsed[[1]][i,3];
		startLine=as.integer(parsed[[1]][i,4]);
		rleValues <- read.table(textConnection(parsed[[1]][i,5]));
		rleO = rle(0);
		rleO$values=as.numeric(rleValues[seq(2,length(rleValues),2)]);
		rleO$lengths=as.integer(rleValues[seq(1,length(rleValues),2)]);
		diffs=inverse.rle(rleO);
		mapping_=c(startLine,startLine+cumsum(diffs[-1]));
		
		basename <- tools::file_path_sans_ext(rnwF);		
		syncF = paste0(basename,".synctex");
		
		compressed <- FALSE
		if (file.exists(syncF)) {
			sf=file(syncF);
		} else{
			syncF <- paste(syncF, ".gz", sep = "")
			if (file.exists(syncF)) {
				compressed <- TRUE
				sf <- gzfile(syncF);
			}
		}
		lines <- try(readLines(syncF, warn = FALSE), silent = TRUE)
		if (inherits(lines, "try-error")) 
			stop(paste(f, "cannot be read, no patching done."))
		close(sf)
		postemble =grep("^Postamble:",lines,perl=T)
		re=paste0("^Input:([^:]+):(.*", basename(texF), ")");		
		toRepl=grep(re,lines[seq(1,postemble)],perl=T);
		inputs=str_match(lines[toRepl],re);
		if (length(inputs)>0){
			tag=inputs[,2];
			inputs[,3]=paste0(tools::file_path_sans_ext(inputs[,3]), ".Rnw");
			lines[toRepl]=paste0("Input:",inputs[,2],":",inputs[,3]);
			re=paste0("^([xkgvh\\$\\(\\[]", tag ,"\\,)(\\d+)([\\,:].*)");
			needReplacment=grep(re,lines[seq(1,postemble)],perl=T);
			toRepl=str_match(lines[needReplacment],re);
			toRepl[,3]=as.character(mapping_[as.integer(toRepl[,3])]);
			newlines=paste0(toRepl[,2],toRepl[,3],toRepl[,4]);
			lines[needReplacment]=newlines;
			if (!compressed) {
				sf=file(syncF,"wb");
			} else{
				sf <- gzfile(syncF,"wb");
			}
			writeLines(lines, sf, sep = "\n")
			close(sf);
			print(paste(length(needReplacment), "patches made to",syncF));
		}else{
			print(paste("No patches made to",syncF));
		}
	}
}

#plot function f with x from xl to xh using either lines 'l' or points 'p'
#f must have single parameter x, can be called with a multi-variable function e.g. profit, like this
#plotf(function(x) {return (profit(x,0,10,5))},0,6,'l'), parameters 0,10,5 remain fixed while x varies
fn$plotf<-function(f,xl,xh,p){
	steps<-300
	xvals<-seq(xl,xh,(xh-xl)/steps)
	yvals<-numeric(steps+1)
	for(i in 1:steps+1)
		yvals[i]<-f(xvals[i])
	plot(xvals,yvals,p)
}
#as plotf but with axis lables
fn$plotfext<-function(f,xl,xh,p,xlab,ylab){
	steps<-300
	xvals<-seq(xl,xh,(xh-xl)/steps)
	yvals<-numeric(steps+1)
	for(i in 1:steps+1)
		yvals[i]<-f(xvals[i])
	plot(xvals,yvals,p,xlab=xlab,ylab=ylab)
}


#Shannon information/ensemble entropy
#call with a single value/vector of probabilities see Mackay (2003, p44)
fn$ent<-function(c){
 result<-0
 for(p in c){
	if(p!=0)
		result<-result+(p*log2(1/p))
	}
 return(result)
}

#c1 P(x1|y1), P(x1)
#c2 P(x2|y2), P(x2)
fn$dist<-function(c1, c2){
	result<-c1[1]*log2(c1[1]/c1[2])
	result<-result+c2[1]*log2(c2[1]/c2[2])
	return(result)
}

#from https://www.r-bloggers.com/update-all-user-installed-r-packages-again/
fn$myInstallAll<-function(){
    install.packages( 
    lib  = lib <- .libPaths()[1],
    pkgs = as.data.frame(installed.packages(lib), stringsAsFactors=FALSE)$Package,
    type = 'source')
}
#binomial function n trials, k successes, p(success)
fn$bin<-function(n,k,p){
	choose(n,k)*p^k*(1-p)^(n-k)
}

#number, decimal places
fn$fndp<-function(n,d=3){
	formatC(n,format='f', digits=d,drop0trailing=TRUE)
}

attach(fn)
