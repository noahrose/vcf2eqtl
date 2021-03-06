plotImbalance <-
function(output,curr,titles=T,usecol=NULL,geneName=NULL,pops=NULL){
	if(is.null(pops)){pops=output$pops}
	if(is.null(pops)) {
		cols=rep('black',ncol(output$genos))
	} else{
		cols=rainbow(length(levels(factor(pops))),v=0.8)[as.numeric(factor(pops))]
	}
	currsnp<-output$res[curr,]
	if(!is.null(usecol)) cols= usecol
	names(cols)<-colnames(output$genos)
	obs<-output$alleleObs[[curr]]
	AIs<-cbind(obs[,2]-obs[,1],obs[,1])
	c1=cols[rownames(obs)]
	if(nrow(AIs)==0){stop('Cannot plot imbalance without ref and alt read depths, check AO and RO VCF fields')}
	plot(1,type='n',xlim=c(0.8,2.2),ylim=c(min(AIs),max(AIs)),axes=F,xlab='Allele',ylab='Read counts')
	axis(1,at=1:2,labels=c(currsnp['REF'],currsnp['ALT']))
	axis(2)
	sapply(1:nrow(AIs),function(ind) lines(AIs[ind,],type='b',lwd=2,col=c1[ind],pch=19))
	if(titles) title('Allelic Imbalace')
	box()
}
