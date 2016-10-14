\name{plotImbalance}
\alias{plotImbalance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotImbalance(output, curr, titles = T, usecol = NULL, geneName = NULL, pops = NULL)
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
  \item{geneName}{
%%     ~~Describe \code{geneName} here~~
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
function (output, curr, titles = T, usecol = NULL, geneName = NULL, 
    pops = NULL) 
{
    if (is.null(pops)) {
        pops = output$pops
    }
    if (is.null(pops)) {
        cols = rep("black", ncol(output$genos))
    }
    else {
        cols = rainbow(length(levels(factor(pops))), v = 0.8)[as.numeric(factor(pops))]
    }
    currsnp <- output$res[curr, ]
    if (!is.null(usecol)) 
        cols = usecol
    legend = F
    AIs <- cbind(output$AImat[curr, output$AIdat == "ref"], output$AImat[curr, 
        output$AIdat == "alt"])
    c1 = cols[!is.na(AIs[, 1])]
    if (length(which(!is.na(AIs[, 1]))) == 0) {
        stop("Cannot plot imbalance without ref and alt read depths, check AO and RO VCF fields")
    }
    AIs <- na.omit(AIs)
    plot(1, type = "n", xlim = c(0.8, 2.2), ylim = c(min(AIs), 
        max(AIs)), axes = F, xlab = "Allele", ylab = "Counts")
    axis(1, at = 1:2, labels = c(currsnp["REF"], currsnp["ALT"]))
    axis(2)
    sapply(1:nrow(AIs), function(ind) lines(AIs[ind, ], type = "b", 
        lwd = 2, col = c1[ind], pch = 19))
    if (titles) 
        title("Allelic Imbalace")
    box()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
