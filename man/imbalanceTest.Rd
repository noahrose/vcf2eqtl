\name{imbalanceTest}
\alias{imbalanceTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
imbalanceTest(AImat, AIdat, fitType = c("parametric", "local", "mean"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AImat}{
%%     ~~Describe \code{AImat} here~~
}
  \item{AIdat}{
%%     ~~Describe \code{AIdat} here~~
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
function (AImat, AIdat, fitType = c("parametric", "local", "mean")) 
{
    cat("reordering into groups with same number of hets...\n")
    counts = AImat
    predictor = do.call(rbind, rep(list(AIdat), nrow(counts)))
    predictor[is.na(counts)] <- NA
    counts[is.na(predictor)] <- NA
    orders <- t(apply(predictor, 1, order))
    c.ord <- t(sapply(1:nrow(orders), function(i) counts[i, orders[i, 
        ]]))
    rownames(c.ord) <- rownames(counts)
    p.ord <- t(sapply(1:nrow(orders), function(i) predictor[i, 
        orders[i, ]]))
    rownames(p.ord) <- rownames(counts)
    notAllNA <- apply(p.ord, 1, function(v) any(!is.na(v)))
    c.ord <- c.ord[notAllNA, ]
    p.ord <- p.ord[notAllNA, ]
    bins <- apply(p.ord, 1, function(v) length(which(!is.na(v))))/2
    bintab <- table(bins)[table(bins) > 1]
    cat("running DESeq...\n")
    res <- do.call(rbind, sapply(names(bintab), DEbin, bins = bins, 
        c.ord = c.ord, p.ord = p.ord, fitType = fitType))
    res <- as.data.frame(res)[rownames(counts), ]
    rownames(res) <- rownames(counts)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
