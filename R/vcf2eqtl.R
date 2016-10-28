vcf2eqtl <-
function(vcf,
expr,
pops=NULL,
minHet=3,
mc.cores=1,
alpha=0.05,
calculateFst=T,
outliers=T,
testDE=F,
all3=T,
hweFilter=T,
hweAlpha=0.05,
covariates=NULL,
propExplained=T,
withinPop=T,
transcripts=NULL,
keepSamples=NULL){
	
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
	AOs<-apply(genoInfo$AO,2,as.numeric)
	ROs<-apply(genoInfo$RO,2,as.numeric)
	AOs[genos!=1]<-NA
	ROs[genos!=1]<-NA
	AImat<-cbind(ROs,AOs)
	colnames(AImat)<-NULL
	rownames(AImat)<-rownames(genos)
	AIdat<-factor(c(rep('ref',ncol(genos)),rep('alt',ncol(genos))),levels=c('ref','alt'))

	cat('organizing SNP info...\n')
	CHROM=as.vector(seqnames(rowRanges(currvcf)))
	POS=as.numeric(as.character(start(ranges(rowRanges(currvcf)))))
	REF=as.character(mcols(rowRanges(currvcf))[,'REF'])
	ALT=as.character(unlist(mcols(rowRanges(currvcf))[,'ALT']))
	AF=unlist(info(currvcf)$AF)
	snpInfo<-cbind(CHROM,POS,REF,ALT,AF)
	rownames(snpInfo)<-rownames(genos)

	cat('organizing expression data...\n')
	#organize and normalize expression data
	if(is.null(transcripts)) transcripts=CHROM
	expr<-expr[rownames(expr)%in%transcripts,]
	voomExpr<-voom(expr)
	rownames(voomExpr$weights)<-rownames(expr)
	currexpr<-as.matrix(voomExpr$E[transcripts,])
	currweights<-as.matrix(voomExpr$weights[transcripts,])
	rownames(currexpr)<-rownames(genos)
	rownames(currweights)<-rownames(genos)
	
	#subset samples if desired
	if(!is.null(keepSamples)){
		cat('subsetting to only keep specified samples and assuming pops correspond to post-filtered samples...\n')
		keep=colnames(genos)%in%keepSamples
		expr<-expr[,keep]
		currexpr<-currexpr[,keep]
		currweights<-currweights[,keep]
		AImat<-AImat[,rep(keep,2)]
		AIdat<-AIdat[rep(keep,2)]
		genos<-genos[,keep]
	}

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
	if(!is.null(minHet)){
		cat('filtering out sites with fewer than',minHet,'heterozygotes...\n')
		numHet<-apply(genos,1,function(v) length(which(v==1)))
		genos<-genos[numHet>=minHet,]
	}

	cat('filtering out sites without allele observations\n')
	AImat<-AImat[rownames(genos),]
	imbalanceInfo<-apply(AImat,1,function(v) any(!is.na(v)))
	genos<-genos[rownames(AImat)[imbalanceInfo],]
	if(nrow(genos)==0){
		stop('No sites left after filtering, check to make sure you have a freebayes VCF of biallelic SNPs with allele observations in it')
	}	
	#subset other data sets after filtering genos
	AImat<-AImat[rownames(genos),]
	currexpr<-currexpr[rownames(genos),]
	snpInfo<-snpInfo[rownames(genos),]
	cat(paste(nrow(genos),'sites left after filtering, testing for eQTL status...\n'))

	cat('allelic imbalance test...\n')
	imb.out<-imbalanceTest(AImat=AImat,AIdat=AIdat)
	
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

	#collect results and calculate p values using Stouffer's method
	res<-cbind(imb.out[,c('stat','pvalue')],2^imb.out$log2FoldChange,assoc.out)
	colnames(res)<-c('AIz','AIp','AIfc','ASSOCz','ASSOCp','ASSOCfc')	
	rownames(res)<-rownames(genos)
	res<-as.data.frame(res)
	res$z<-(res[,'AIz']+res[,'ASSOCz'])/(2**.5)
	# res$z[is.na(res$AIz)&!is.na(res$ASSOCz)]<-res$ASSOCz[is.na(res$AIz)&!is.na(res$ASSOCz)]
	# res$z[!is.na(res$AIz)&is.na(res$ASSOCz)]<-res$ASSOCz[!is.na(res$AIz)&is.na(res$ASSOCz)]
	res$p<-sapply(res$z,function(val) min(pnorm(val,lower.tail=T),pnorm(val,lower.tail=F))*2)
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
		wc.out<-MakeDiploidFSTMat(t(genos),rownames(genos),pops)
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
	
	return(list(res=res,snpContigExpr=currexpr,genos=genos,AImat=AImat,AIdat=AIdat,globalFst=globalFst,pops=pops,DEres= DEres))
}
