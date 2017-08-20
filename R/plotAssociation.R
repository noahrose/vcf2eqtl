plotAssociation <-
function(output,curr,titles=T,plotLegend=T,legendPos=NULL,snpLab=NULL,usecol=NULL,geneName=NULL,lty=2,lwd=2,pops=NULL,residual=T){
	if(is.null(pops)){pops=output$pops}
	if(is.null(pops)) {
		plotLegend=F
		cols=rep('black',ncol(output$genos))
	}
	ex=output$snpContigExpr[curr,]
	if(residual){
		ex<-scale(lm(ex~pops)$residuals)
	}
	gt=output$genos[curr,]
	currsnp<-output$res[curr,]
	if(!is.null(pops)) cols=rainbow(length(levels((pops))),v=0.8)[as.numeric((pops))]
	if(!is.null(usecol)) {cols=usecol; plotLegend=F}
	xlab=paste(currsnp$CHROM,'pos.', currsnp['POS'])
	if(!is.null(snpLab)) xlab=snpLab
	plot(ex~jitter(gt,0.5),col=cols,pch=19,axes=F, ylab='Log2 CPM', xlab=xlab,xlim=c(-.1,2.1))
	abline(lm(ex~gt),lwd=lwd,lty=lty)
	axis(2)
	legpos='topleft'
	if(currsnp$ASSOCz<0) legpos='topright'
	if(!is.null(legendPos)) legpos=legendPos
	if(plotLegend) legend(legpos,fill=rainbow(length(levels(pops)),v=0.8),legend=levels(pops),bty='n')
	axis(1,at=c(0,1,2),labels=c(paste(currsnp['REF'],currsnp['REF'],sep=''),
		paste(currsnp['REF'],currsnp['ALT'],sep=''),
		paste(currsnp['ALT'],currsnp['ALT'],sep='')))
	if(titles) title('eQTL Association')
	box()
	if(!is.null(geneName)) mtext(geneName,1,2)
}
