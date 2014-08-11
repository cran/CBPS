\name{CBPS}
\alias{CBPS}
\alias{CBPS.fit}
\alias{print.CBPS}
\title{Covariate Balancing Propensity Score (CBPS) Estimation}
\description{
  \code{CBPS} estimates propensity scores such that both covariate balance and prediction of treatment assignment are maximized.   The method, therefore, avoids an iterative process between model fitting and balance checking and implements both simultaneously.
For cross-sectional data, the method can take continuous treatments and treatments with a control (baseline) condition and either 1, 2, or 3 distinct treatment conditions.}
\usage{
	  CBPS(formula, data, na.action, ATT = NULL, method = "over",
	       iterations = NULL, standardize = TRUE, twostep = TRUE, ...)
	  CBPS.fit(treat, X, ATT, X.bal = X, method, iterations, 
			   standardize, twostep, ...)
}
\arguments{
  \item{formula}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
  \item{na.action}{A function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset.}
  \item{ATT}{Set to \code{TRUE} to find the average treatment effect on the treated,
	         and \code{FALSE} to find the average treatment effect.  Default is \code{TRUE}.  For 3- and 4-valued treatments as well as continuous treatments , only the ATE is available.}
  \item{treat}{A vector of treatment assignments.}
  \item{X}{A covariate matrix.}
  \item{X.bal}{A matrix to be balanced.}
  \item{method}{Choose "over" to fit an over-identified model that combines the propensity score and covariate balancing conditions; choose "exact" to fit a model that only contains the covariate balancing conditions.}
  \item{iterations}{An optional parameter for the maximum number of iterations for the optimization.  Default is 1000.}
  \item{standardize}{Default is \code{TRUE}, which gets inverse propensity score weights (IPW) as described in Hirano, Imbens, and Ridder (2003).   Set to \code{FALSE} to return Horvitz-Thompson weights.  Does not affect weights for continuous treatments.}
  \item{twostep}{Default is \code{TRUE} for a two-step estimator, which will run substantially faster than continuous-updating.  Set to \code{FALSE} to use the continuous-updating estimator described by Imai and Ratkovic (2014).}
  \item{...}{Other parameters to be passed through to \code{optim()}}
}
\details{Fits covariate balancing propensity scores.}
\value{
  \item{coefficients}{A named vector of coefficients}
  \item{sigmasq}{The sigma-squared value, for continuous treatments only.}
  \item{residuals}{The working residuals from the final iteration of CBPS}
  \item{fitted.values}{The fitted propensity score.}
  \item{rank}{The numeric rank of the fitted model}
  \item{family}{The family object used.}
  \item{deviance}{Minus twice the log-likelihood of the CBPS fit.  Will be lower than the maximum likelihood deviance.}
  \item{weights}{The optimal weights.  For the marginal structural model and MultiBins, returns a list in which stabilized ("stabilized") weights, unstabilized ("unstabilized") weights, and unconditional treatment probabilities ("unconditional") are all available.  Weights are Horvitz-Thompson if standardize is set to \code{FALSE}, and IPW if standardize is set to \code{TRUE}.}
  \item{y}{The treatment vector used}
  \item{x}{The covariate matrix}
  \item{model}{The model frame}
  \item{converged}{Convergence value.  Returned from the call to \code{optim()}.}
  \item{call}{The matched call}
  \item{formula}{The formula supplied.}
  \item{data}{The data argument.}
  \item{J}{The J-statistic at convergence.}
  \item{mle.J}{The J-statistic for the parameters from maximum likelihood estimation.}
  \item{bal}{The balance loss at convergence.}
  \item{mle.bal}{The balance loss for the paramaters from maximum likelihood estimation.}
  \item{df}{The degrees of freedom.}
  \item{var}{The covariance matrix, evaluated numerically from \code{optim()}.}
}
\references{Imai, Kosuke and Marc Ratkovic.  2014. ``Covariate Balancing Propensity Score.'' Journal of the Royal Statistical Society, Series B (Statistical Methodology). \url{http://imai.princeton.edu/research/CBPS.html}				
			Fong, Christian and Kosuke Imai.  ``Covariate Balancing Propensity Scores for General Treatment Regimes.''  Unpublished Manuscript. \url{http://imai.princeton.edu/research/files/CBGPS.pdf}
			}
\author{Christian Fong, Marc Ratkovic, and Kosuke Imai; The CBPS function is based on the code for version 2.15.0 of the glm function implemented in the stats package, originally written by Simon Davies.  This documenation is likewise modeled on the documentation for glm and borrows its language where the arguments and values are the same.}


\seealso{\link{summary.CBPS}}

\examples{
\dontrun{
###
### Example: propensity score matching
###

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS via logistic regression
fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), 
			data = LaLonde, ATT = TRUE)
summary(fit)
## matching via MatchIt: one to one nearest neighbor with replacement
library(MatchIt)
m.out <- matchit(treat ~ fitted(fit), method = "nearest", data = LaLonde, 
                 replace = TRUE)


### Example: propensity score weighting 
###
## Simulation from Kang and Shafer (2007).
set.seed(123456)
n <- 500
X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 0.25*X[,3] + 0.1 * X[,4]))
treat <- rbinom(n, 1, prop)
y <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)

##Estimate CBPS with a misspecificied model
X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, (X[,1]*X[,3]/25+.6)^3, 
               (X[,2]+X[,4]+20)^2)
fit1 <- CBPS(treat ~ X.mis, ATT = FALSE)
summary(fit1)
	
## Horwitz-Thompson estimate
mean(treat*y/fit1$fitted.values)
## Inverse propensity score weighting
sum(treat*y/fit1$fitted.values)/sum(treat/fit1$fitted.values)

rm(list=c("y","X","prop","treat","n","X.mis","fit1"))

### Example: Covariate Balancing Propensity Score - 3 Treatments
##
set.seed(1)
n<-1000
k<-3
X<-mvrnorm(n,rep(0,k), Sigma=matrix(c(5,.5,-.03,.5,1,-.27,-.03,-.27,1),3,3))
beta<-matrix(2*k,k,2)
prob<-cbind((1+exp(X\%*\%beta[,1])+exp(X\%*\%beta[,2]))^-1,
			 exp(X\%*\%beta[,1])*(1+exp(X\%*\%beta[,1])+exp(X\%*\%beta[,2]))^-1,
			 exp(X\%*\%beta[,2])*(1+exp(X\%*\%beta[,1])+exp(X\%*\%beta[,2]))^-1)
treat.latent<-runif(n)
treat<-ifelse(treat.latent < prob[,1], 1, 
			  ifelse(treat.latent < (prob[,1] + prob[,2]), 2, 3))
			  
gamma<-rnorm(ncol(X))
y<-X\%*\%gamma + ifelse(treat == 2, 2, ifelse(treat==3, -1, 0)) + rnorm(1000)

cbpsfit3<-CBPS(treat ~ X)
outfit3<-lm(y ~ as.factor(treat) + X, weights = cbpsfit3$weights)

### Example: Continuous Treatment
set.seed(1)
n<-1000
X<-cbind(mvrnorm(n, mu = rep(0,4), Sigma = diag(4)), rbinom(n,1,runif(1)), rbinom(n,10,0.7))
beta<-rnorm(ncol(X)+1, sd=5)
noise<-rbinom(1,10,0.5)*.50
treat<-cbind(1,X)\%*\%beta+rnorm(n, sd = noise)
prob<-dnorm(treat, mean = cbind(1,X)\%*\%beta, sd = noise, log = TRUE)

treat.effect<-rnorm(1)
effectbeta<-rnorm(ncol(X))
y<-rbinom(n, 1, (1+exp(mean(treat.effect*treat + X\%*\%effectbeta) - 
							treat.effect*treat - X\%*\%effectbeta))^-1)

cbpsfitc<-CBPS(treat ~ X)
outfitc<-lm(y ~ treat + X, weights = cbpsfitc$weights) 
}
}