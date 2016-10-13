propExplain <-
function(curr,currexpr,genos,pops){
	ex=currexpr[curr,]
	gt=genos[curr,]
	aov1<-summary(aov(ex~pops))[[1]]	
	aov2<-summary(aov(ex~gt+pops))[[1]]	
	popDiffExplained=max(1-aov2['pops','Sum Sq']/aov1['pops','Sum Sq'],0)
	return(popDiffExplained)
}
