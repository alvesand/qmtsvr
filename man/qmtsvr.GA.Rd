\name{qmtsvr.GA}
\alias{qmtsvr.GA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QMTSVR hyperparameters optimization using a Genetic Algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
```{r setup}
pheno = data_ex$mydata
X = data_ex$QTL

hyper_st = list(c("C",0.1,2,64), c("eps",0.001,0.1,64), c("b1",0.1,7,64))

hyper_mt = list(c("C",0.1,2,64), c("eps",0.001,0.1,64),
                c("b1",0.1,7,64), c("b2",0.1,7,64), c("b3",0.1,7,64),
                c("b4",0.1,7,64), c("b5",0.1,7,64), c("b6",0.1,7,64),
                c("r1",-0.9,0.9,64),c("r2",-0.9,0.9,64),c("r3",-0.9,0.9,64),
                c("r4",-0.9,0.9,64),c("r5",-0.9,0.9,64), c("r6",-0.9,0.9,64))


sol_svr = qmtsvr.GA(Y = pheno[100:300,5:7], X = X[100:300,], hyper = hyper_mt, vartype = "continuous", dopar = F,
                    ngen = 30, popsize = 40, mut_rate = 0.04, cross_rate = 0.95,
                    elitism = 2, cost = "cor", tsize = 5, val_pop = "cross", nfolds = 3,lambda=0,
                    vardiag=T,MRCR = "dynamic", verbose = T)


sol_svr$GA_parameters
plot.GA(sol_svr)

```

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
