\name{DEbin}
\alias{DEbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DEbin(bin, bins, c.ord, p.ord, fitType = c("parametric", "local", "mean"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bin}{
%%     ~~Describe \code{bin} here~~
}
  \item{bins}{
%%     ~~Describe \code{bins} here~~
}
  \item{c.ord}{
%%     ~~Describe \code{c.ord} here~~
}
  \item{p.ord}{
%%     ~~Describe \code{p.ord} here~~
}
  \item{fitType}{
%%     ~~Describe \code{fitType} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bin, bins, c.ord, p.ord, fitType = c("parametric", 
    "local", "mean")) 
{
    numBin = length(which(bins == bin))
    currc <- c.ord[bins == bin, ]
    currp <- p.ord[bins == bin, ][1, ]
    if (numBin < 20) 
        fitType = "mean"
    nafilt = !is.na(currp)
    currc <- currc[, nafilt]
    currp <- currp[nafilt]
    form = as.formula(~currp)
    currdat = data.frame(currp = factor(currp))
    cat(paste("testing sites with", table(currdat$currp)[1], 
        "heterozygotes...\n"))
    cds <- DESeqDataSetFromMatrix(currc, currdat, form)
    cds <- DESeq(cds, fitType = fitType)
    res <- results(cds)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
