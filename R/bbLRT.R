bbLRT<-function(dat){
	full.out<-suppressWarnings(optim(c(1,1),bbfull,data=dat))
	null.out<-suppressWarnings(optim(c(1,1),bbnull,data=dat))
	p.out<-pchisq(2*(null.out$value-full.out$value),df=1,lower.tail=F)
	mu<-full.out$par[1]/sum(full.out$par)
	l2fc<-log2(mu/(1-mu))
	z=sign(0.5-mu)*qnorm(p)
	return(c(z=z,p.value=p.out, l2fc=l2fc))
}