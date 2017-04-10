bbLRT<-function(dat){
	full.out<-suppressWarnings(optim(c(1,1),bbFull,data=dat))
	null.out<-suppressWarnings(optim(c(1,1),bbNull,data=dat))
	p.out<-pchisq(2*(null.out$value-full.out$value),df=1,lower.tail=F)
	mu<-full.out$par[1]/sum(full.out$par)
	l2fc<-log2(mu/(1-mu))
	return(c(mu=mu,p.value=p.out, l2fc=l2fc))
}