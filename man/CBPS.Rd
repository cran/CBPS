\name{CBPS}
\alias{CBPS}
\alias{CBPS.fit}
\alias{print.CBPS}
\title{Covariate Balancing Propensity Score (CBPS) Estimation}
\description{
  \code{CBPS} and \code{CBPS.fit} estimate propensity score such that both covariate balance and prediction of treatment assignment are maximized.  The method, therefore, avoids an iterative process between model fitting and balance checking and implements both simultaneously.
}
\usage{
	  CBPS(formula, data, na.action, ATT, k = 0, ...)
	  CBPS.fit(treat, X, ATT, X.bal = X, k, ...)
}
\arguments{
  \item{formula}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
    \item{na.action}{A function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset.}
  \item{ATT}{Set to \code{TRUE} to find the average treatment effect on the treated,
	         and \code{FALSE} to find the average treatment effect.  Default is \code{TRUE}.}
  \item{treat}{A vector of treatment assignments.}
  \item{X}{A covariate matrix.}
  \item{X.bal}{A matrix to be balanced. Currently not accepted.}
  \item{k}{The model dimensionality. 0 is the default, and is estimated from X. Useful with ill-conditioned matrices.}
  \item{...}{Other parameters to be passed through to \code{optim()}}
}
\details{Fits covariate balancing propensity scores.}
\value{
  \item{coefficients}{A named vector of coefficients}
  \item{residuals}{The working residuals from the final iteration of CBPS}
  \item{fitted.values}{The fitted mean values, obtained by transforming the linear predictors by the inverse of the binomial function}
  \item{rank}{The numeric rank of the fitted model}
  \item{family}{The family object used.}
  \item{deviance}{Minus twice the log-likelihood of the CBPS fit.  Will be lower than the maximum likelihood deviance.}
  \item{weights}{The CBPS balancing weights calculated from the propensity scores.}
  \item{y}{The treatment vector used}
  \item{x}{The covariate matrix}
  \item{model}{The model frame}
  \item{converged}{Convergence value.  Returned from the call to \code{optim()}.}
  \item{call}{The matched call}
  \item{formula}{The formula supplied.}
  \item{data}{The data argument.}
  \item{J}{The J-statistic at convergence.}
  \item{df}{The degrees of freedom.}
  \item{var}{The covariance matrix, evaluated numerically from \code{optim()}.}
}
\references{Imai, Kosuke and Marc Ratkovic.  ``Covariate Balancing Propensity Score.'' Unpublished Manuscript, Princeton University. \url{http://imai.princeton.edu/research/CBPS.html}}

\author{Marc Ratkovic, Kosuke Imai, and Christian Fong; The CBPS function is based on the code for version 2.15.0 of the glm function implemented in the stats package, originally written by Simon Davies.  This documenation is likewise modeled on the documentation for glm and borrows its language where the arguments and values are the same.}

\seealso{\link{summary.CBPS}}

\examples{
###
### Example: propensity score matching
###

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS via logistic regression
fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), data = LaLonde, ATT = TRUE)
summary(fit)
## matching via MatchIt: one to one nearest neighbor with replacement
library(MatchIt)
m.out <- matchit(treat ~ 1, distance = fitted(fit), method = "nearest", data = LaLonde, replace = TRUE)


### Example: propensity score weighting 
###
## Simulation from Kang and Shafer (2007).
set.seed(123456)
n <- 1000
X <- mvrnorm(n, mu = rep(0, 4), Sigma = diag(4))
prop <- 1 / (1 + exp(X[,1] - 0.5 * X[,2] + 0.25*X[,3] + 0.1 * X[,4]))
treat <- rbinom(n, 1, prop)
y <- 210 + 27.4*X[,1] + 13.7*X[,2] + 13.7*X[,3] + 13.7*X[,4] + rnorm(n)

##Estimate CBPS with a misspecificied model
X.mis <- cbind(exp(X[,1]/2), X[,2]*(1+exp(X[,1]))^(-1)+10, (X[,1]*X[,3]/25+.6)^3, (X[,2]+X[,4]+20)^2)
fit1 <- CBPS(treat ~ X.mis, ATT = TRUE)
	
## Horwitz-Thompson estimate
mean(treat*y/fit1$fitted.values)
## Inverse probability weighting
sum(treat*y/fit1$fitted.values)/sum(treat/fit1$fitted.values)

}