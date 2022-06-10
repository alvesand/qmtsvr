\name{qmtsvr.dist}
\alias{qmtsvr.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Euclidean Distances
%%  ~~ Fast computation of an N x N (Weighted) Euclidean Matrix ~~
}
\description{
%%  ~~ This function uses a numeric n x p Matrix X (e.g. marker matrix)
for computing an n x n Euclidean Distance Matrix (EDM). The function can also accomodate different weighting factors for the variables in X.~~
}
\usage{
qmtsvr.dist(x=X, w=NULL, u=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{x    matrix to compute the EDM
w     a vector of weights for computing the Euclidean distances (optional)
u     the power of the final EDM values (default = 2)
  \item{x}{
%%     ~~Describe \code{x} here~~
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
## Create toy data
values = rnorm(5000, 0, 1)
X = matrix(values, 1000, 50)
D = qmtsvr.dist (x=X)

dim(D)
head(D[,1:10])

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
