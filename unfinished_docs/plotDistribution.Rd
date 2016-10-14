\name{plotDistribution}
\alias{plotDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDistribution(output, curr, titles = T, usecol = NULL, plotLegend = T, legendPos = "topright", xlab = "Population", pops = NULL)
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
  \item{usecol}{
%%     ~~Describe \code{usecol} here~~
}
  \item{plotLegend}{
%%     ~~Describe \code{plotLegend} here~~
}
  \item{legendPos}{
%%     ~~Describe \code{legendPos} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
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
function (output, curr, titles = T, usecol = NULL, plotLegend = T, 
    legendPos = "topright", xlab = "Population", pops = NULL) 
{
    if (is.null(pops)) {
        if (!is.null(output$pops)) {
            pops = output$pops
        }
        else {
            stop("Need populations to display allele distributions between populations...")
        }
    }
    cols = rainbow(length(levels(pops)), v = 0.8)[1:length(levels(pops))]
    if (!is.null(usecol)) 
        cols = usecol
    legend = F
    agg <- aggregate(output$genos[curr, ] ~ pops, FUN = mean)
    xlim = c(1, length(levels(pops)))
    xlim = mean(xlim) + (xlim - mean(xlim)) * 1.25
    ylim = c(min(agg[, 2]), max(agg[, 2]))/2
    ylim = mean(ylim) + (ylim - mean(ylim)) * 1.25
    plot(agg[, 2]/2, col = cols, pch = 19, cex = 2, xlab = xlab, 
        ylab = "Alternate allele freq.", axes = F, xlim = xlim, 
        ylim = ylim)
    axis(1, at = 1:length(levels(pops)), levels(pops))
    axis(2)
    lines(agg[, 2]/2, lty = 2, lwd = 2)
    box()
    if (titles) 
        title("Allele distribution")
    if (plotLegend) {
        legend(legendPos, fill = cols, legend = levels(pops), 
            bty = "n")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
