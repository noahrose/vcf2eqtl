makeHierfstat <-
function(mat,pops,remove.names=F){
	newmat=mat
	newmat[mat==0]<-11
	newmat[mat==1]<-12
	newmat[mat==2]<-22
	newmat<-as.data.frame(cbind(factor(pops),t(newmat)))
	if(remove.names) {rownames(newmat)=NULL;colnames(newmat)=c('pops',1:(ncol(newmat)-1))}
	return(newmat)
}
