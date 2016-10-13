DEbin <-
function(bin,bins,c.ord,p.ord,fitType=c('parametric','local','mean')){
	numBin= length(which(bins==bin))
	currc<-c.ord[bins==bin,]
	currp<-p.ord[bins==bin,][1,]
	if(numBin<20) fitType='mean'
	nafilt=!is.na(currp)
	currc<-currc[,nafilt]
	currp<-currp[nafilt]
	form=as.formula(~currp)
	currdat=data.frame(currp=factor(currp))
	cat(paste('testing sites with',table(currdat$currp)[1],'heterozygotes...\n'))	
	cds<-DESeqDataSetFromMatrix(currc,currdat,form)
	cds<-DESeq(cds,fitType=fitType)
	res<-results(cds)
	return(res)
}
