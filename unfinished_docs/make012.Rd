\name{make012}
\alias{make012}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make012(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
%%     ~~Describe \code{mat} here~~
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
function (mat) 
{
    row_names = rownames(mat)
    col_names = colnames(mat)
    allowed = c("0|0", "0/0", "0/1", "1/0", "0|1", "1|0", "1/1", 
        "1|1")
    if (any(!mat \%in\% allowed)) {
        stop(paste("unexpected genotype(s):", paste(mat[!mat \%in\% 
            allowed], collapse = ","), "-- are you sure these are biallelic SNPs?"))
    }
    newmat <- mat
    newmat[mat == "0/0" | mat == "0|0"] <- 0
    newmat[mat == "0/1" | mat == "1/0" | mat == "0|1" | mat == 
        "1|0"] <- 1
    newmat[mat == "1/1" | mat == "1|1"] <- 2
    newmat <- apply(newmat, 2, as.numeric)
    rownames(newmat) <- row_names
    colnames(newmat) <- col_names
    return(newmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
