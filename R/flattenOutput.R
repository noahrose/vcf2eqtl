flattenOutput <-
function(output){
	ex<-output$snpContigExpr
	colnames(ex)<-paste(colnames(ex),'_expr',sep='')
	gt<-output$genos
	colnames(gt)<-paste(colnames(gt),'_geno',sep='')
	return(cbind(output$res,ex,gt))
}
