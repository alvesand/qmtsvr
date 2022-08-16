\name{qmtsvr.dist}
\alias{qmtsvr.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Weighted Euclidean Distances
%%  ~~ Fast computation of an N x N (Weighted) Euclidean Matrix ~~
}
\description{
%%  ~~ This function uses a numeric n x p Matrix X (e.g. marker matrix)
This function is used for computing an \emph{n x n} Euclidean Distance Matrix (EDM). It can also accomodate different weighting factors for the variables in X.
}
\usage{
qmtsvr.dist(x=X, w=NULL, u=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{A \emph{n x p} matrix}
\item{w}{A vector of size \emph{p} containing the weighting factors to be assigned to the features in \bold{X} during the EDM computation (\emph{Optional}).}
\item{u}{the power of the final EDM values (\emph{default =} 2).}
\item{verbose}{A logical argument indicating if the time elapsed for computing the EDM must be printed (\emph{default = F})}
\item{vardiag}{A logical argument specifying if EDM diagonal elements must vary (\emph{default=}\code{F}). If \code{vardiag}=\code{T}, the zero values in the diagonal of the EDM will be replaced with: \deqn{1-(\frac{d}{\sqrt d * \sqrt m})} in which \bold{d} is the diagonal of \eqn{XX^T} and \emph{m} is the median of \bold{d}.}
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
Anderson Alves (carvalhoalve@wisc.edu)

webpage: \url{http://alvesand.netlify.app/}
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
