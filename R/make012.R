make012 <-
function(mat){
	row_names=rownames(mat)
	col_names=colnames(mat)
	allowed=c('0|0','0/0','0/1','1/0','0|1','1|0','1/1','1|1','.')
	if(any(!mat%in%allowed)){
		stop(paste('unexpected genotype(s):',paste(mat[!mat%in%allowed],collapse=','),'-- are you sure these are biallelic SNPs?'))
	}
	newmat<-mat
	newmat[mat=='0/0'|mat=='0|0']<-0
	newmat[mat=='0/1'|mat=='1/0'|mat=='0|1'|mat=='1|0']<-1
	newmat[mat=='1/1'|mat=='1|1']<-2
	newmat[mat=='.']<-NA
	newmat<-apply(newmat,2,as.numeric)
	rownames(newmat)<-row_names
	colnames(newmat)<-col_names
	return(newmat)
}
