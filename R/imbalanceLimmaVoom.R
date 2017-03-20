imbalanceLimmaVoom <-
function(currsnp,AImat,AIdat,AIwt){
	res<-summary(lm(AImat[currsnp,]~AIdat,weights=AIwt[currsnp,]))$coefficients[2,c(3,4,1)]
	return(res)
}
