bbNull<-function(par,data){
	-sum(dbetabinom.ab(data$x,data$size,par[1],par[1],log=T))
}