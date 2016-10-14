\name{vcf2eqtl}
\alias{vcf2eqtl}
\alias{associationTest}
\alias{DEbin}
\alias{flattenOutput}
\alias{imbalanceTest}
\alias{make012}
\alias{makeHierfstat}
\alias{plotAssociation}
\alias{plotCline}
\alias{plotImbalance}
\alias{plotQTL}
\alias{propExplan}
\alias{vcf2eqtl-package}

\title{
Function to detect eQTLs using vcf and expression matrix from a transcriptome reference
}
\description{
A simple pipeline that takes a freebayes vcf file of biallelic SNPs and an 
expression matrix in DESeq format and uses two different tests to detect eQTLs.
First, we test for allelic imbalance in heterozygous individuals using DESeq2.
Overdispersion in read counts is modeled by testing all sites with a given number of
heterozygotes together in a single DESeq2 run. Next, we use limma/voom to log transform
and calculate precision weights for expression values and use lm to test for associations
between genotype and expression. P values are calculated using Stouffer's method.
Fst and Fst outliers (OutFLANK), Differential expression between populations (DESeq2), 
and the proportion of population difference in expression explained by eQTLs can be 
calculated.
}
\usage{
vcf2eqtl(vcf, expr, pops = NULL, minHet = 3, mc.cores = 1, alpha = 0.05, calculateFst = T, testDE = F, all3 = T, hweFilter = T, hweAlpha = 0.05, covariates = NULL, propExplained = T, withinPop = T, transcripts = NULL, keepSamples = NULL)
}
\arguments{
  \item{vcf}{
A freebayes-generated vcf of biallelic SNPs from RNA-Seq reads mapped to a transcriptome reference
}
  \item{expr}{
An expression matrix, same format as DESeq. Rows are contigs, Columns are samples (same order as in vcf), rownames are contig names
}
  \item{pops}{
A factor vector of which population each sample belongs to
}
  \item{minHet}{
Minimum number of heterozygotes to test for allelic imbalance
}
  \item{mc.cores}{
Number of cores to use for multithreading
}
  \item{alpha}{
FDR for calling eQTLs and differential expression
}
  \item{calculateFst}{
Boolean; do you want to calculate Fst?
}
  \item{testDE}{
Boolean; test for differential expression between pops? Uses DeSeq2 LRT
}
  \item{all3}{
Boolean; require all three genotypes to test for eQTL at a site?
}
  \item{hweFilter}{
Boolean; filter out sites out of Hardy-Weinberg equilibrium?
}
  \item{hweAlpha}{
HWE P value in at least one population to filter out site
}
  \item{covariates}{
Matrix of covariates to include in association test
}
  \item{propExplained}{
Boolean; calculate proportion of population difference in expression explained by eQTL?
}
  \item{withinPop}{
Boolean; do association tests within populations to account for population structure?
}
  \item{transcripts}{
Boolean; list of transcripts that each site falls into; useful if you mapped to a genome.
However, if you mapped to a genome, better tools may exist (e.g. WASP) to test for eQTLs.
}
  \item{keepSamples}{
Vector of sample names to retain for analysis, assumes that your pops vector corresponds
to post-filtered samples. Does not change sample order for retained samples.
}
}
\details{
This function is useful for finding eQTLs among biallelic SNPs in RNA-Seq data mapped to
a transcriptome reference, and is most useful for less-resouced non-model organisms.
If you have a good genome and phased data, more powerful methods like WASP will likely
outperform this function.
}

res = res, snpContigExpr = currexpr, genos = genos, 
        AImat = AImat, AIdat = AIdat, globalFst = globalFst, 
        pops = pops, DEres = DEres

\value{
  \item{res} data frame of SNPs with eQTL status and other info (as specified) like Fst etc
  \item{snpContigExpr} voom-transformed expression for contig associated with each test SNP
  \item{genos} 012 matrix of genotypes
  \item{AImat} matrix of allele observations for ref and alt within heterozygotes
  \item{AIdat} vector of which columns in AImat are ref and alt
  \item{globalFst} globalFst estimate, if calculated
  \item{pops} population vector
  \item{DEres} results from DESeq2, if calculated
}
\references{
Genomic responses to environmental variation and change in the ocean. PhD Dissertation,
Noah Rose, Stanford University.
}
\author{
Noah Rose
}
\note{
}


\seealso{
}
\examples{
expr <- read.delim('expr.txt',row.names=1)
vcf <- 'biallelic_snps.vcf'
metadata<-read.delim('meta.txt')
pops=metadata$pops
vcf2eqtl(vcf,expr,pops)
}
\keyword{eQTL}
\keyword{VCF}
