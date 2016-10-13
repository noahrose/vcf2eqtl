plotDistribution <-
function(output,curr,titles=T,usecol=NULL,plotLegend=T,legendPos='topright',xlab='Population',pops=NULL){
	if(is.null(pops)) {
		if(!is.null(output$pops)){
			pops=output$pops
		} else{
		stop('Need populations to display allele distributions between populations...')
		}
	}
	cols=rainbow(length(levels(pops)),v=0.8)[1:length(levels(pops))]
	if(!is.null(usecol)) cols= usecol; legend=F
	agg<-aggregate(output$genos[curr,]~pops,FUN=mean)
	xlim=c(1,length(levels(pops)))
	xlim=mean(xlim)+(xlim-mean(xlim))*1.25
	ylim=c(min(agg[,2]),max(agg[,2]))/2
	ylim=mean(ylim)+(ylim-mean(ylim))*1.25
	plot(agg[,2]/2,col=cols,pch=19,cex=2,xlab=xlab,ylab='Alternate allele freq.',axes=F,xlim=xlim,ylim=ylim)
	axis(1,at=1:length(levels(pops)),levels(pops))
	axis(2)
	lines(agg[,2]/2,lty=2,lwd=2)
	box()
	if(titles) title('Allele distribution')
	if(plotLegend) {legend(legendPos,fill=cols,legend=levels(pops),bty='n')}
}
