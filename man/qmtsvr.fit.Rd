\name{qmtsvr.fit}
\alias{qmtsvr.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting (\emph{Quasi}) Multi-trait Support Vector Regression methods
%%  ~~ Provides an API for fitting QMTSVR models using the kernlab package~~
}
\description{
This function fits a QMTSVR model with with user pre-defined hyper parameters (requires the package kernlab). Returns an \emph{n x t} YHAT matrix.

%%  ~~ Provides an API for fitting QMTSVR models using the kernlab package ... ~~
}
\usage{
qmtsvr.fit = function(Y=Y, X=X, set_hyper=set_hyper, w=NULL, D=NULL, verbose=F, vardiag=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Y}{an \emph{n x t} numeric matrix of correlated response variables, where \emph{N} is the number of observations and \emph{t} is the number of traits (missing values are coded as \code{NA}).}
\item{X}{A \emph{N x p} matrix of features (\emph{e.g.,} SNP markers ), in which \emph{p} is the number of features.}
\item{set_hyper}{a list specifying the QMTSVR model hyperparameters}
\item{w}{A list containing the trait-specific and trait-common weighting factors to be assigned to the features in \bold{X} during the EDM computation (\emph{Optional}).}
\item{D}{a list with all pre-computed EDM. It can be passed as an argument alternatively to the matrix \bold{X} (\emph{Optional}).}
\item{verbose}{A logical argument indicating if the function log must be printed at the console.}
\item{vardiag}{A logical argument specifying if diagonal elements of the Euclidean Distance Matrix (EDM) must vary (\emph{default=F}). If \code{vardiag}=\code{T}, the zero values in the diagonal of the EDM will be replaced with:   \deqn{1-(\frac{d}{\sqrt d * \sqrt m})} in which \bold{d} is the diagonal of \eqn{XX^T} and \emph{m} is the median of \bold{d}.This quantity yields a variation around 1 for the diagonal elements of the RBF kernel, computed as \eqn{K = e^{-b*D^{2}}}, where \emph{b} is a hyperparameter constant and \bold{D} is the Euclidean Distance Matrix. The diagonal elements in \bold{K} correlate with those from the \bold{G} matrix presented in VanRaden, 2008.}

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

webpage:\url{http://alvesand.netlify.app/}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[=dest]{qmtsvr.dist},\link[=dest]{qmtsvr.GA}}
}
\examples{

library(qmtsvr)

#Load dataset
pheno = data_ex$mydata[1:300,5:7] #For illustration, only first 300 observations are used
X = data_ex$SNP[1:300,] # Load SNP markers
obs = pheno #copy the observations matrix
pheno[200:300,1] = NA #omit the last 100 observations

#Set arbitrary hyperparameters values for a multi-task SVR model with 3 traits
hyper_list = list("C" = 0.5, "eps" = 0.03, "b1" = 4.5, "b2"= 3.67, "b3" = 6.73,"b4" = 1.25,  "b5" = 0.95, "b6" = 1.89, "r1" = -0.40,
            "r2" = 0.53, "r3" = 0.32)

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

#Predict the values
YHAT = qmtsvr.fit(Y=pheno, X = X, set_hyper = hyper_list, verbose = T, vardiag = T)

#Correlate the predictions with the observed values
cor(YHAT[200:300,1],obs[200:300,1])
#compute Root-mean squared error
mean(sqrt((YHAT[200:300,1] - obs[200:300,1])^2))

##Optimize the hyperparameters using the qmtsvr.GA function
##Define hyperparameters range for the model
hyper_mt = list(c("C",0.1,5,128), c("eps",0.001,0.3,128),
                c("b1",0.1,7,128), c("b2",0.1,7,128), c("b3",0.1,7,128),
                c("b4",0.1,7,128), c("b5",0.1,7,128), c("b6",0.1,7,128),
                c("r1",-0.9,0.9,128),c("r2",-0.9,0.9,128),c("r3",-0.9,0.9,128))

#Run the Genetic algorithm with 3-fold cross-validation
mt_svr = qmtsvr.GA(Y = pheno, X = X[1:300,], hyper = hyper_mt, ngen = 30, popsize = 50, mut_rate = 0.02, cross_rate = 0.95,
                    elitism = 2, cost = "cor", tsize = 5, val_pop = "cross", nfolds = 3, vardiag=T, verbose = T)

#Use the hyperparameters list obtained from mt_svr to predict the unobserved values
YHAT2 = qmtsvr.fit(Y=pheno, X = X, set_hyper = mt_svr$set_hyper, verbose = T, vardiag = T)

#Correlate the predictions with the observed values
cor(YHAT2[200:300,1],obs[200:300,1])
#compute Root-mean squared error
mean(sqrt((YHAT2[200:300,1] - obs[200:300,1])^2))

#plot the results
par(mfrow=c(1,2))
plot(YHAT[200:300,1],obs[200:300,1], pch = 19, lty = 2, xlab = "Predicted values", ylab = "Observed values", col = rgb(0.7,0.5,0.20,0.7), main = "Arbitrary hyperparameters")



plot(YHAT2[200:300,1],obs[200:300,1], pch = 19, lty = 2, xlab = "Predicted values",ylab = "Observed values", col = rgb(0.1,0.3,0.40,0.5), main = "Optimized hyperparameters")


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
