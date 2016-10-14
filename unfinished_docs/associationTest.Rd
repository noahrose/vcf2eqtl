\name{associationTest}
\alias{associationTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
associationTest(curr, currexpr, currweights, genos, covariates = NULL, withinPop = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{curr}{
%%     ~~Describe \code{curr} here~~
}
  \item{currexpr}{
%%     ~~Describe \code{currexpr} here~~
}
  \item{currweights}{
%%     ~~Describe \code{currweights} here~~
}
  \item{genos}{
%%     ~~Describe \code{genos} here~~
}
  \item{covariates}{
%%     ~~Describe \code{covariates} here~~
}
  \item{withinPop}{
%%     ~~Describe \code{withinPop} here~~
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
function (curr, currexpr, currweights, genos, covariates = NULL, 
    withinPop = T) 
{
    gt = genos[curr, ]
    components = "gt"
    if (withinPop) 
        components = c("pops", "gt")
    if (!is.null(covariates)) 
        components = c("covariates", components)
    form <- as.formula(paste("ex~", paste(components, collapse = "+")))
    ex = currexpr[curr, ]
    if (is.null(currweights)) {
        lm.out <- summary(lm(form))$coefficients
    }
    else {
        wt = currweights[curr, ]
        lm.out <- summary(lm(form, weights = wt))$coefficients
    }
    if (!"gt" \%in\% rownames(lm.out)) {
        warning(paste("could not fit", curr, " -- skipping and returning NA..."))
        return(rep(NA, 3))
    }
    lm.res <- lm.out["gt", c("t value", "Pr(>|t|)")]
    fc = 2^(lm.out["gt", "Estimate"] * 2)
    return(c(lm.res, fc))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
