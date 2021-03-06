# vcf2eqtl

##A pipeline for eQTL detection from a transcriptome reference
- The pipeline uses a simple counts matrix and VCF input file
- Both of these file types are standard parts of RNA-Seq experiments

##Install (R 3.3.1 or greater):
```R
source("https://bioconductor.org/biocLite.R")
biocLite(c('devtools','DESeq2','limma','VariantAnnotation','qvalue'))
library(devtools)
install_github('whitlock/OutFLANK')
install_github('noahrose/vcf2eqtl')
```

###Call variants with freebayes and filter for biallelic SNPs
```
freebayes -f ref.fa *.bam > fb.vcf
cat fb.vcf \
| vcfnulldotslashdot \
| grep -vF './.' \
| vcffilter -f "TYPE = snp & QUAL > 30 & AF > 0.1 & AF < 0.9 & NUMALT = 1" -g "DP > 10" \
> biallelic_snps.vcf
```

###You will also need a gene expression counts matrix of the same format used for DESeq
- Many ways to do this, including:
  * samtools idxstats
  * htseq-count
  * Trinity/RSEM
  * StringTie
  * etc...

###Then load into R

```
expr <- read.delim('expr.txt',row.names=1)
vcf <- 'biallelic_snps.vcf'
metadata<-read.delim('meta.txt')
pops=metadata$pops
vcf2eqtl(vcf,expr,pops)
```