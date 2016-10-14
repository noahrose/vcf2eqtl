\name{plotAssociation}
\alias{plotAssociation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotAssociation(output, curr, titles = T, plotLegend = T, legendPos = NULL, snpLab = NULL, usecol = NULL, geneName = NULL, lty = 2, lwd = 2, pops = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
%%     ~~Describe \code{output} here~~
}
  \item{curr}{
%%     ~~Describe \code{curr} here~~
}
  \item{titles}{
%%     ~~Describe \code{titles} here~~
}
  \item{plotLegend}{
%%     ~~Describe \code{plotLegend} here~~
}
  \item{legendPos}{
%%     ~~Describe \code{legendPos} here~~
}
  \item{snpLab}{
%%     ~~Describe \code{snpLab} here~~
}
  \item{usecol}{
%%     ~~Describe \code{usecol} here~~
}
  \item{geneName}{
%%     ~~Describe \code{geneName} here~~
}
  \item{lty}{
%%     ~~Describe \code{lty} here~~
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
}
  \item{pops}{
%%     ~~Describe \code{pops} here~~
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
function (output, curr, titles = T, plotLegend = T, legendPos = NULL, 
    snpLab = NULL, usecol = NULL, geneName = NULL, lty = 2, lwd = 2, 
    pops = NULL) 
{
    if (is.null(pops)) {
        pops = output$pops
    }
    if (is.null(pops)) {
        plotLegend = F
        cols = rep("black", ncol(output$genos))
    }
    ex = output$snpContigExpr[curr, ]
    gt = output$genos[curr, ]
    currsnp <- output$res[curr, ]
    if (!is.null(pops)) 
        cols = rainbow(length(levels((pops))), v = 0.8)[as.numeric((pops))]
    if (!is.null(usecol)) {
        cols = usecol
        plotLegend = F
    }
    xlab = paste(currsnp$CHROM, "pos.", currsnp["POS"])
    if (!is.null(snpLab)) 
        xlab = snpLab
    plot(ex ~ jitter(gt, 0.5), col = cols, pch = 19, axes = F, 
        ylab = "Log2 Normalized Counts", xlab = xlab, xlim = c(-0.1, 
            2.1))
    abline(lm(ex ~ gt), lwd = lwd, lty = lty)
    axis(2)
    legpos = "topleft"
    if (currsnp$z < 0) 
        legpos = "topright"
    if (!is.null(legendPos)) 
        legpos = legendPos
    if (plotLegend) 
        legend(legpos, fill = rainbow(length(levels(pops)), v = 0.8), 
            legend = levels(pops), bty = "n")
    axis(1, at = c(0, 1, 2), labels = c(paste(currsnp["REF"], 
        currsnp["REF"], sep = ""), paste(currsnp["REF"], currsnp["ALT"], 
        sep = ""), paste(currsnp["ALT"], currsnp["ALT"], sep = "")))
    if (titles) 
        title("eQTL Association")
    box()
    if (!is.null(geneName)) 
        mtext(geneName, 1, 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
