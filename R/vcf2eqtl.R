vcf2eqtl <-
function(vcf,
expr,
pops=NULL,
minHet=3,
minHetDP=10,
mc.cores=1,
alpha=0.05,
calculateFst=T,
outliers=T,
testDE=F,
all3=F,
hweFilter=T,
hweAlpha=0.05,
covariates=NULL,
propExplained=T,
withinPop=T,
format='bcftools',
transcripts=NULL){
	
	if(is.null(pops)){
		calculateFst=F
		testDE=F
		propExplained=F
		withinPop=F
	}
	#extract SNP genotype and SNP depth data from vcf
	globalFst=NULL
	cat('reading vcf...\n')
	currvcf<-suppressWarnings(readVcf(vcf,genome='curr'))
	genoInfo<-geno(currvcf)
	genos<-make012(genoInfo$GT)
	
	if(format=='freebayes'){
		AOs<-apply(genoInfo$AO,2,as.numeric)
		ROs<-apply(genoInfo$RO,2,as.numeric)
	} else if (format=='bcftools'){
		ad<-apply(geno(currvcf)$AD,c(1,2),function(l) l[[1]])
		ROs<-ad[1,,]
		AOs<-ad[2,,]
	} else{
		stop('format must be either bcftools or freebayes')
	}
	
		BOs<-AOs+ROs
		AOs[genos!=1]<-NA
		BOs[genos!=1]<-NA
		AOs[BOs<minHetDP]<-NA
		BOs[BOs<minHetDP]<-NA
	

	rownames(AOs)<-rownames(genos)
	rownames(BOs)<-rownames(genos)
	cat('getting reference and alternate observations...\n')
	alleleObs<-list()
	for(i in 1:nrow(AOs)){
		alleleObs[[rownames(genos)[i]]]<-na.omit(data.frame(row.names=colnames(genos),x=AOs[i,],size=BOs[i,]))
	}	

	cat('organizing SNP info...\n')
	CHROM=as.vector(seqnames(rowRanges(currvcf)))
	POS=as.numeric(as.character(start(ranges(rowRanges(currvcf)))))
	REF=as.character(mcols(rowRanges(currvcf))[,'REF'])
	ALT=as.character(unlist(mcols(rowRanges(currvcf))[,'ALT']))
	AF=as.numeric(unlist(info(currvcf)$AC))/as.numeric(unlist(info(currvcf)$AN))
	snpInfo<-cbind(CHROM,POS,REF,ALT,AF)
	rownames(snpInfo)<-rownames(genos)

	cat('organizing expression data...\n')
	#organize and normalize expression data
	if(is.null(transcripts)) transcripts=CHROM
	expr<-expr[rownames(expr)%in%transcripts,]
	designMatrix<-NULL
	dmatComponents<-NULL
	if(withinPop) dmatComponents ='pops'
	if(!is.null(covariates)) dmatComponents =c('covariates', dmatComponents)
	if(!is.null(dmatComponents)){
		form<-as.formula(paste('~',paste(dmatComponents,collapse='+')))
		designMatrix<-model.matrix(form)	
	}
	voomExpr<-voom(expr,design=designMatrix)
	rownames(voomExpr$weights)<-rownames(expr)
	currexpr<-as.matrix(voomExpr$E[transcripts,])
	currweights<-as.matrix(voomExpr$weights[transcripts,])
	rownames(currexpr)<-rownames(genos)
	rownames(currweights)<-rownames(genos)

	#if desired, filter genotypes	
	if(all3){
		cat('filtering for sites with all three genotypes...\n')
		num_gts<-apply(genos,1,function(v) length(table(v)))
		genos<-genos[num_gts==3,]
	}
	if(hweFilter){
		cat('filtering out sites out of HWE in at least one population...\n')
		hwepops=pops
		hwe=rep(1,nrow(genos))
		if(is.null(pops)) hwepops=rep(1,ncol(genos))
		for(pop in unique(hwepops)){
			hwe<-pmin(hwe,apply(genos[,which(pops==pop)],1,function(v) HWExact(table(factor(v,levels=c(0,1,2))),verbose=F)$pval))
		}
		genos<-genos[hwe>hweAlpha,]
	}

	cat('filtering out sites without minimum number of informative heterozygotes\n')
	alleleObs<-alleleObs[rownames(genos)]
	imbalanceInfo<-unlist(lapply(alleleObs,function(df) nrow(df)>minHet))
	alleleObs<-alleleObs[imbalanceInfo]
	genos<-genos[names(alleleObs),]
	if(nrow(genos)==0){
		stop('No sites left after filtering, check to make sure you have a freebayes VCF of biallelic SNPs with allele observations in it')
	}	
	#subset other data sets after filtering genos
	currexpr<-currexpr[rownames(genos),]
	snpInfo<-snpInfo[rownames(genos),]
	cat(paste(nrow(genos),'sites left after filtering, testing for eQTL status...\n'))

	#single-threaded test
	if(mc.cores==1){	
		cat('imbalance test...\n')
		imb.out<-do.call(rbind,lapply(alleleObs,bbLRT))
	} else{
	#multithreaded test
		cat('mulithreaded imbalance test...\n')
		imb.out<-do.call(rbind,mclapply(alleleObs,bbLRT,mc.cores=mc.cores))
	}
	
	#single-threaded test
	if(mc.cores==1){	
		cat('association test...\n')
		assoc.out<-t(sapply(rownames(genos),associationTest,
		currexpr=currexpr,currweights=currweights,genos=genos,withinPop=withinPop))
	} else{
	#multithreaded test
		cat('mulithreaded association test...\n')
		assoc.out<-do.call(rbind,mclapply(rownames(genos),associationTest,
			currexpr=currexpr,currweights=currweights,genos=genos,
			withinPop=withinPop,mc.cores=mc.cores))
	}

	#collect results and calculate p values using Fisher's method
	res<-cbind(imb.out,assoc.out)
	colnames(res)<-c('AImu','AIp','AIlog2fc','ASSOCz','ASSOCp','ASSOClog2fc')	
	rownames(res)<-rownames(genos)
	res<-as.data.frame(res)
	#compare alternate homozygotes for fc
	combineP<-function(v){
		if(NA%in%v) return(NA)
		return(sumlog(v)$p)
	}
	res$ASSOClog2fc<-2*res$ASSOClog2fc
	res$p<-apply(cbind(res$AIp,res$ASSOCp),1,combineP)
	res$padj<-p.adjust(res$p,method='BH')
	res$AIpadj<-p.adjust(res$AIp,method='BH')
	res$ASSOCpadj<-p.adjust(res$ASSOCp,method='BH')
	res$eQTL<-res$padj<alpha
	res<-cbind(snpInfo,res)
	res$REF<-as.character(res$REF)
	res$ALT<-as.character(res$ALT)
	res$POS<-as.numeric(as.character(res$POS))
	
	#calculate Fst and call outliers using OutFLANK
	if(calculateFst){
		cat('calculating Fst...\n')
		genosOutFLANK=genos
		genosOutFLANK[is.na(genosOutFLANK)]<-9
		wc.out<-MakeDiploidFSTMat(t(genosOutFLANK),rownames(genos),pops)
		res$Fst=wc.out$FST
		res$FstNum<-wc.out$T1
		res$FstDen<-wc.out$T2
		res$He<-wc.out$He
		globalFst=sum(res$FstNum)/sum(res$FstDen)
		if(outliers){
			fl.out<-OutFLANK(wc.out,NumberOfSamples=ncol(genos),qthreshold=alpha)$results
			res$FstOutlier=fl.out$OutlierFlag
			res$FstOutlierP=fl.out$pvaluesRightTail
		}
	}
	
	#test for differential expression using DESeq2
	DEres =NULL
	if(testDE){
		cat('testing for differential expression with DESeq2...\n')
		currdat<-data.frame(pops=pops)
		cds<-DESeqDataSetFromMatrix(expr,currdat,~pops)
		cds<-DESeq(cds,test='LRT',reduced=~1)
		DEres <-results(cds)
		res$DE<-DEres[res$CHROM,'padj']<alpha
		res$DEp<-DEres[res$CHROM,'pvalue']
		res$DEpadj<-DEres[res$CHROM,'padj']
	}
	
	#calculate reduction in population differentiation after accounting for eQTL
	if(propExplained){
		cat('calculating proportion of population differences explained by eQTLs...\n')
		if(mc.cores==1){
			propE<-sapply(rownames(genos)[which(res$eQTL)],propExplain,currexpr=currexpr,genos=genos,pops=pops)
		} else{
		cat('\tmulithreading...\n')
			propE<-unlist(mclapply(rownames(genos)[which(res$eQTL)],propExplain,
				currexpr=currexpr,genos=genos,pops=pops,mc.cores=mc.cores))
		}
		res$popDiffExplained<-NA
		res$popDiffExplained[which(res$eQTL)]=propE
	}
	
	return(list(res=res,snpContigExpr=currexpr,genos=genos,alleleObs=alleleObs,globalFst=globalFst,pops=pops,DEres= DEres))
}
