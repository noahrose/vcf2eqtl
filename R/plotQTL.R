plotQTL <-
function(output,curr,
	titles=T,
	plotLegend=T,
	ourpar=T,
	legendPos='topright',
	snpLab=NULL,
	col=NULL,
	geneName=NULL,
	lty=2,
	lwd=2,
	cline=NULL,
	pops=NULL,
	plotDistribution=T){
	if(ourpar) {
		par(tck=-0.01,mgp=c(1.2,0.2,0),mar=c(3,3,2,1),bty='l')
		if(plotDistribution) {
			par(mfrow=c(1,3))	
		} else{
			par(mfrow=c(1,2))
		}
	}
	plotAssociation(output,curr,titles,plotLegend=F,legendPos,snpLab,col,geneName,lty,lwd,pops=pops)
	plotImbalance(output,curr,titles,col,geneName,pops=pops)
	if(plotDistribution){
		if(is.null(cline)){
			plotDistribution(output,curr,titles,col,plotLegend,legendPos,pops=pops)
		} else{
			plotCline(cline,output,curr,titles,col,plotLegend,legendPos)
		}
	}
}
