\name{qmtsvr.GA}
\alias{qmtsvr.GA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{QMTSVR hyperparameters optimization using a Genetic Algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
qmtsvr.GA(Y = Y, X = X, w = NULL, D = NULL, hyper = hyper, ngen = ngen, popsize = popsize, mut_rate = mut_rate, cross_rate = cross_rate, cost = "cor", tsize = tsize, val_pop = "cross", nfolds = 3)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Y}{A \emph{N x t} matrix of response variables, where N is the number of observations and t is the number of traits (missing values are coded as NA).}
\item{X}{A \emph{N x p} matrix of features (\emph{e.g.,} SNP markers ), in which \emph{p} is the number of features.}
\item{w}{A list containing the trait-specific and trait-common weighting factors to be assigned to the features in \bold{X} during the EDM computation (\emph{Optional}).}
\item{D}{a list with all pre-computed EDM. It can be passed as an argument alternatively to the matrix \bold{X} (\emph{Optional}).}
\item{hyper}{A \emph{list} object containing the name specification and range  for each hyperparameter to be optimized in the QMTSVR model}
\item{ngen}{The number of generations to run in the Genetic Algorithm.}
\item{popsize}{The population size (PS) to grow, it indicates the number of models tested per generation.}
\item{mut_rate}{The mutation rate (\emph{e.g.,} 0.05) of the genetic algorithm. It is the probability of a single bit (or gene) on the binary array changing randomly, implying slight modifications to the current model}
\item{cross_rate}{The crossing over (CR) rate for the genetic algorithm. The CR controls the rate that a child's chromosome will result from the crossing-over of two parents instead of being an identical copy of one of them.}
\item{tsize}{The tournament size to be used in the GA(\emph{default=}5).For each child chromosome to be created, \emph{tsize} individuals are drawn at random from the current population and the one with the highest fs is selected for integrating the pair of crossing chromosomes }
\item{elitism}{The number of models with the best performance in the current generation to keep for the next generarion (\emph{default=}1).}
\item{cost}{The cost function to be optimized in the GA. The values accepted are: \code{"cor"},\code{"rmse"} and \code{"mae"} for Perason correlation between observed and predicted values, \emph{Root-mean squared error} and \emph{Mean absolute error}, respectively. The current version only accepts cost functions for continuous target variables.}
\item{w}{A vector of size \emph{p} with the weighting factors to be assigned to the features in \bold{X} for the EDM computation (Optional).}
\item{tgTrait}{A constant for indexing the column where the target variable is located at Y (\emph{default=}1).}
\item{val_pop}{A character object indicating the validation strategy adopted in the GA. \code{"cross"} indicates cross-validation, \code{"custom"} indicates that the cost for an specific group must be optimized and \code{"closest"} indicates that the target population to be optimized will be composed by the closest individuals (based on the Euclidean distance of X) to those with NA values in the target variable.}
\item{nfolds}{If val_pop = \code{"cross"}, \emph{nfolds} is the number of folds to be used in the cross-validation.}
\item{k}{If val_pop = \code{"closest"}, \emph{k} is the number closest neighbors to be identified in the EDM for each NA value.}
\item{custom_val}{If val_pop = \code{"custom"}, \emph{custom_val} is a vector of size \emph{N} with 0 (train) and 1 (validation), indicating which individual will compose the validation population in the GA.}
\item{vartype}{The target variable class. The current version only accepts \code{"continuous"}.}
\item{vardiag}{A logical argument specifying if diagonal elements of the Euclidean Distance Matrix (EDM) must vary (\emph{default=F}). If \code{vardiag}=\code{T}, the zero values in the diagonal of the EDM will be replaced with:   \deqn{1-(\frac{d}{\sqrt d * \sqrt m})} in which \bold{d} is the diagonal of \eqn{XX^T} and \emph{m} is the median of \bold{d}.This quantity yields a variation around 1 for the diagonal elements of the RBF kernel, computed as \eqn{K = e^{-b*D^{2}}}, where \emph{b} is a hyperparameter constant and \bold{D} is the Euclidean Distance Matrix. The diagonal elements in \bold{K} correlate with those from the \bold{G} matrix presented in VanRaden, 2008.}
\item{MRCR}{\code{"fixed"}(\emph{default}) or \code{"dynamic"}, if the second option is chosen, the GA will gradatively reduce CR and increase MR when the algorithm get stucked in local minima.}
\item{verbose}{A logical argument indicating if the GA output must be printed at the console.}
\item{dopar}{A logical argument indicating if the GA will be performed using a parallelized process. It will speed up the GA runtime. \bold{Attention: in the current version it will only work in Linux machines}.}
\item{ncores}{If \emph{dopar}=\code{T}, the number of cores to be used must be specified (\emph{default=}number of cores - 1).}
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

library(qmtsvr)

#Load dataset
pheno = data_ex$mydata
X = data_ex$SNP

##Define hyperparameters range for a single-task SVR
hyper_st = list(c("C",0.01,5,128), c("eps",0.001,0.2,128), c("b1",0.1,7,128))
##The regularization parameter ("C") will range from 0.01 to 5, with 128 values tested within this interval
##The epsilon constant ("eps") will range from 0.001 to 0.2, with 128 values tested within this interval
##The bandwidth parameter ("b1") will range from 0.1 to 7, with 128 values tested within this interval

##Run the GA for a single-task SVR (as illustration only the first 200 values are used)
## The validation strategy used is a 3-fold cross-validation
st_svr = qmtsvr.GA(Y = as.matrix(pheno[1:200,5]), X = X[1:200,], hyper = hyper_st, ngen = 30, popsize = 100, mut_rate = 0.05, cross_rate = 0.95, elitism = 2, cost = "mae", tsize = 5, val_pop = "cross", nfolds = 3, vardiag=T, verbose = T)

#Get the list with optimized hyperparameters
st_svr$set_hyper

##Define hyperparameters range for a multi-task SVR with 3 variables
hyper_mt = list(c("C",0.1,5,128), c("eps",0.001,0.3,128),
                c("b1",0.1,7,128), c("b2",0.1,7,128), c("b3",0.1,7,128),
                c("b4",0.1,7,128), c("b5",0.1,7,128), c("b6",0.1,7,128),
                c("r1",-0.9,0.9,128),c("r2",-0.9,0.9,128),c("r3",-0.9,0.9,128))
#The sequence of bandwith parameters (b) according to the traits (t1, t2, t3) are:
#     t1 t2 t3
# t1: b1 b2 b3
# t2:  . b4 b5
# t3:  .  . b6

#The sequence of weights r according to the traits (t1, t2, t3) are:
#     t1 t2 t3
# t1:  . r1 r2
# t2:  .  . r3
# t3:  .  .  .

##Run the GA for a multi-task SVR with 3 response variables
mt_svr = qmtsvr.GA(Y = pheno[1:200,5:7], X = X[1:200,], hyper = hyper_mt, ngen = 30, popsize = 100, mut_rate = 0.02, cross_rate = 0.95,
                   elitism = 2, cost = "cor", tsize = 5, val_pop = "cross", nfolds = 3,vardiag=T, verbose = T)


#Get the list with optimized hyperparameter values
mt_svr$set_hyper

##Plot GA history
par(mfrow=c(1,2))
plot.GA(st_svr)
plot.GA(mt_svr)

##Run the same MTSVR model but using the last 50 individuals as the GA validation population
custom_val = c(rep(0,150),rep(1,50))

mt_svr2 = qmtsvr.GA(Y = pheno[1:200,5:7], X = X[1:200,], hyper = hyper_mt, ngen = 30, popsize = 100, mut_rate = 0.02, cross_rate = 0.95,
                   elitism = 2, cost = "cor", tsize = 5, val_pop = "custom",
                   custom_val = custom_val, vardiag=T, verbose = T)

plot.GA(mt_svr2)

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
