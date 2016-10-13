imbalanceTest <-
function(AImat,AIdat,fitType=c('parametric','local','mean')){
	cat('reordering into groups with same number of hets...\n')
	counts=AImat
	predictor=do.call(rbind,rep(list(AIdat),nrow(counts)))
	predictor[is.na(counts)]<-NA
	counts[is.na(predictor)]<-NA
	orders<-t(apply(predictor,1,order))
	c.ord<-t(sapply(1:nrow(orders),function(i) counts[i,orders[i,]]))
	rownames(c.ord)<-rownames(counts)
	p.ord<-t(sapply(1:nrow(orders),function(i) predictor[i,orders[i,]]))
	rownames(p.ord)<-rownames(counts)
	notAllNA<-apply(p.ord,1,function(v) any(!is.na(v)))
	c.ord<-c.ord[notAllNA,]
	p.ord<-p.ord[notAllNA,]
	bins<-apply(p.ord,1,function(v) length(which(!is.na(v))))/2
	bintab<-table(bins)[table(bins)>1]
	cat('running DESeq...\n')
	res<-do.call(rbind,sapply(names(bintab),DEbin,bins=bins,c.ord=c.ord,p.ord=p.ord,fitType=fitType))
	res<-as.data.frame(res)[rownames(counts),]
	rownames(res)<-rownames(counts)
	return(res)
}
