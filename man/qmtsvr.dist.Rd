\name{qmtsvr.dist}
\alias{qmtsvr.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Euclidean Distances
%%  ~~ Fast computation of an N x N (Weighted) Euclidean Matrix ~~
}
\description{
%%  ~~ This function uses a numeric n x p Matrix X (e.g. marker matrix)
Function used for computing an n x n Euclidean Distance Matrix (EDM). The function can also accomodate different weighting factors for the variables in X.
}
\usage{
qmtsvr.dist(x=X, w=NULL, u=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
x :    an n x p matrix

w :     a vector of weights for computing the Euclidean distances (optional)

u :     the power of the final EDM values (default = 2)

verbose :     prints the time elapsed for computing the EDM (default = F)

vardiag :      if vardiag is set True, the zeros in the diagonal of the EDM will be replaced with 1-(d/(sqrt(d)*sqrt(m)), in which d is the diagonal of tcrossprod(x) and m is the median of this diagonal.This quantity yields a variation around 1 for the diagonal elements of the RBF kernel computed as K = exp(-b*EDM), where b is a hyperparameter constant. The diagonal elements in K correlate with the diagonal elements from the G matrix presented in VanRaden, 2008. (default = F)
{
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
values = rnorm(50, 0, 1)
X = matrix(values, 5, 10)
D = qmtsvr.dist (x=X)
D
## Give more weigth for the collumns 1 and 5
w = rep(1,50)
w[1] = 5
w[5] = 15.5
D2 = qmtsvr.dist (x=X, w = w)
D2
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
