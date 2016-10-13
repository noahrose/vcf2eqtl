plotCline <-
function(cline,output,curr,titles=T,usecol=NULL,plotLegend=T,legendPos='topright',xlab='Latitude'){
	clinefac<-factor(cline)
	cols=rainbow(length(levels(clinefac)),v=0.8)[1:length(levels(clinefac))]
	if(!is.null(usecol)) cols= usecol; legend=F
	agg<-aggregate(output$genos[curr,]~cline,FUN=mean)
	xlim=c(min(agg[,1]),max(agg[,1]))
	xlim=mean(xlim)+(xlim-mean(xlim))*1.25
	ylim=c(min(agg[,2]),max(agg[,2]))/2
	ylim=mean(ylim)+(ylim-mean(ylim))*1.25
	plot(agg[,1],agg[,2]/2,col=cols,pch=19,cex=2,xlab=xlab,ylab='Alternate allele freq.',axes=F,xlim=xlim,ylim=ylim)
	axis(1)
	axis(2)
	lines(agg[,1],agg[,2]/2,lty=2,lwd=2)
	box()
	if(titles) title('Allele distribution')
}
