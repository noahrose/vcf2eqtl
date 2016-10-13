associationTest <-
function(curr,currexpr,currweights,genos,covariates=NULL,withinPop=T){
	gt=genos[curr,]
	components='gt'
	if(withinPop) components=c('pops','gt')
	if(!is.null(covariates)) components=c('covariates',components)
	form<-as.formula(paste('ex~',paste(components,collapse='+')))
	ex= currexpr[curr,]
	if(is.null(currweights)) {
		lm.out<-summary(lm(form))$coefficients
	} else{
		wt= currweights[curr,]
		lm.out<-summary(lm(form,weights=wt))$coefficients
	}
	if(!'gt'%in%rownames(lm.out)) {
		warning(paste('could not fit',curr,' -- skipping and returning NA...'))
		return(rep(NA,3))
	}
	lm.res<-lm.out['gt',c('t value','Pr(>|t|)')]
	fc=2^(lm.out['gt','Estimate']*2)
	return(c(lm.res,fc))
}
