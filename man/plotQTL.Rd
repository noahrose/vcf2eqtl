\name{plotQTL}
\alias{plotQTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotQTL(output, curr, titles = T, plotLegend = T, ourpar = T, legendPos = "topright", snpLab = NULL, col = NULL, geneName = NULL, lty = 2, lwd = 2, cline = NULL, pops = NULL, plotDistribution = T)
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
  \item{ourpar}{
%%     ~~Describe \code{ourpar} here~~
}
  \item{legendPos}{
%%     ~~Describe \code{legendPos} here~~
}
  \item{snpLab}{
%%     ~~Describe \code{snpLab} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
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
  \item{cline}{
%%     ~~Describe \code{cline} here~~
}
  \item{pops}{
%%     ~~Describe \code{pops} here~~
}
  \item{plotDistribution}{
%%     ~~Describe \code{plotDistribution} here~~
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
function (output, curr, titles = T, plotLegend = T, ourpar = T, 
    legendPos = "topright", snpLab = NULL, col = NULL, geneName = NULL, 
    lty = 2, lwd = 2, cline = NULL, pops = NULL, plotDistribution = T) 
{
    if (ourpar) {
        par(tck = -0.01, mgp = c(1.2, 0.2, 0), mar = c(3, 3, 
            2, 1), bty = "l")
        if (plotDistribution) {
            par(mfrow = c(1, 3))
        }
        else {
            par(mfrow = c(1, 2))
        }
    }
    plotAssociation(output, curr, titles, plotLegend = F, legendPos, 
        snpLab, col, geneName, lty, lwd, pops = pops)
    plotImbalance(output, curr, titles, col, geneName, pops = pops)
    if (plotDistribution) {
        if (is.null(cline)) {
            plotDistribution(output, curr, titles, col, plotLegend, 
                legendPos, pops = pops)
        }
        else {
            plotCline(cline, output, curr, titles, col, plotLegend, 
                legendPos)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
