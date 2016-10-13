\name{propExplain}
\alias{propExplain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
propExplain(curr, currexpr, genos, pops)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{curr}{
%%     ~~Describe \code{curr} here~~
}
  \item{currexpr}{
%%     ~~Describe \code{currexpr} here~~
}
  \item{genos}{
%%     ~~Describe \code{genos} here~~
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
function (curr, currexpr, genos, pops) 
{
    ex = currexpr[curr, ]
    gt = genos[curr, ]
    aov1 <- summary(aov(ex ~ pops))[[1]]
    aov2 <- summary(aov(ex ~ gt + pops))[[1]]
    popDiffExplained = max(1 - aov2["pops", "Sum Sq"]/aov1["pops", 
        "Sum Sq"], 0)
    return(popDiffExplained)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
