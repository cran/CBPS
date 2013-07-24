\name{CBPS}
\alias{CBPS}
\alias{CBPS.fit}
\alias{CBMSM.fit}
\alias{print.CBPS}
\title{Covariate Balancing Propensity Score (CBPS) Estimation}
\description{
  \code{CBPS} estimates propensity scores such that both covariate balance and prediction of treatment assignment are maximized.   The method, therefore, avoids an iterative process between model fitting and balance checking and implements both simultaneously.
For cross-sectional data, the method can take treatments with a control (baseline) condition and either 1, 2, or 3 distinct treatment conditions.  With longitudinal data, the method returns marginal structural model weights that can be entered directly into a linear model.  The method also handles multiple binary treatments adminstered concurrently.}
\usage{
	  CBPS(formula, data, na.action, ATT = NULL, method = "over",
	       type="propensity", iterations = NULL, standardize = FALSE, ...)
	  CBPS.fit(treat, X, ATT, X.bal = X, method, iterations, standardize, ...)
	  CBMSM.fit(treat, X, time, method, MultiBin.fit = FALSE, standardize, ...)
}
\arguments{
  \item{formula}{For the cross-sectional case (type="propensity"), an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.  For the marginal structural model (type="MSM"), a list of formulas.  The function assumes that there is one formula for each time, and they are ordered from the first time to the last time.}
  \item{data}{An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
    \item{na.action}{A function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of options, and is na.fail if that is unset.}
  \item{ATT}{Set to \code{TRUE} to find the average treatment effect on the treated,
	         and \code{FALSE} to find the average treatment effect.  Default is \code{TRUE}.  For 3 and 4-valued treatments and class "MSM", only the ATE is available.}
  \item{treat}{A vector of treatment assignments.}
  \item{X}{A covariate matrix.}
  \item{X.bal}{A matrix to be balanced. Currently not accepted.}
   \item{method}{"over" for fitting an over-identified model that combines the propensity score and covariate balancing conditions; "exact" for fitting a model that only contains the covariate balancing conditions.}
   \item{type}{Either "propensity" for a single time period and a treatment variable that takes 2, 3, or 4 unique values;"MSM" for a marginal structural model, with multiple time periods; and "MultiBin" for multiple binary treatments at the same time period.}
   \item{time}{A vector giving the time period for each observation; only used for type "MSM"}
   \item{MultiBin.fit}{A parameter for whether the multiple binary treatments occur concurrently (FALSE) or over consecutive time periods (TRUE) as in a marginal structural model.}
   \item{iterations}{An optional parameter for the maximum number of iterations for the optimization.  Default is 1000 for binary and 10000 otherwise.}
   \item{standardize}{Set to \code{TRUE} to get inverse propensity score weights (IPW) as described in Hirano, Imbens, and Ridder (2003).  Default is \code{FALSE}, which returns Horvitz-Thompson weights.}
  \item{...}{Other parameters to be passed through to \code{optim()}}
}
\details{Fits covariate balancing propensity scores.}
\value{
  \item{coefficients}{A named vector of coefficients}
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
  \item{df}{The degrees of freedom.}
  \item{var}{The covariance matrix, evaluated numerically from \code{optim()}.}
}
\references{Imai, Kosuke and Marc Ratkovic.  Forthcoming. ``Covariate Balancing Propensity Score.'' Journal of the Royal Statistical Society, Series B (Statistical Methodology). \url{http://imai.princeton.edu/research/CBPS.html}		
	
	Imai, Kosuke and Marc Ratkovic.  ``Robust Estimation of Inverse Probability Weights for Marginal Structural Models.'' Unpublished Manuscript, Princeton University. \url{http://imai.princeton.edu/research/MSM.html} 
	
	Fong, Christian.  ``Covariate Balancing Propensity Scores for Multi-Valued Treatment Specifications.''  Unpublished Manuscript.
	}

\author{Marc Ratkovic, Kosuke Imai, and Christian Fong; The CBPS function is based on the code for version 2.15.0 of the glm function implemented in the stats package, originally written by Simon Davies.  This documenation is likewise modeled on the documentation for glm and borrows its language where the arguments and values are the same.}

\seealso{\link{summary.CBPS}}

\examples{
\dontrun{
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
## Inverse probability weighting
sum(treat*y/fit1$fitted.values)/sum(treat/fit1$fitted.values)

rm(list=c("y","X","prop","treat","n","X.mis","fit1"))

### Example: Covariate Balancing Marginal Structural Model
###
## Simulation from Imai and Ratkovic (2013).


data1<-MSMdata(200)
attach(data1)
formulas.msm<-list(c(treat.1~X1, treat.2~X2, treat.3~X3))
msm1<-CBPS(formulas.msm,method="over",type="MSM")

#Unweighted model
lm1<-lm(y~treat.1+treat.2+treat.3)
summary(lm1)

#CBMSM weighted model
lm2<-lm(y~treat.1+treat.2+treat.3,weights=msm1$weights)
summary(lm2)

detach(data1)



### Example: Covariate Balancing Propensity Score - 3 Treatments
##

set.seed(1)
library(MASS)
X<-cbind(1,mvrnorm(1000,c(0,0,0), Sigma=matrix(c(5,.5,-.03,.5,1,-.27,-.03,-.27,1),3,3)))
beta<-matrix(rnorm(8),4,2)
prob<-cbind((1+exp(X\%*\%beta[,1])+exp(X\%*\%beta[,2]))^-1,
				 exp(X\%*\%beta[,1])*(1+exp(X\%*\%beta[,1])+exp(X\%*\%beta[,2]))^-1,
				 exp(X\%*\%beta[,2])*(1+exp(X\%*\%beta[,1])+exp(X\%*\%beta[,2]))^-1)
X<-X[,2:4]
treat.latent<-runif(1000)
treat<-ifelse(treat.latent < prob[,1], 1, ifelse(treat.latent < (prob[,1] + prob[,2]), 2, 3))
fit3<-CBPS(treat ~ X, ATT = FALSE)



### Example: Multiple Binary Treatments Administered at the Same Time

n<-200
k<-4
set.seed(1040)
X1<-cbind(1,matrix(rnorm(n*k),ncol=k))

betas.1<-betas.2<-betas.3<-c(2,4,4,-4,3)/5


probs.1<-probs.2<-probs.3<-(1+exp(-X1 \%*\% betas.1))^-1

treat.1<-rbinom(n=length(probs.1),size=1,probs.1)
treat.2<-rbinom(n=length(probs.2),size=1,probs.2)
treat.3<-rbinom(n=length(probs.3),size=1,probs.3)

X1<-X1[,-1]

formulas.multi<-list(c(treat.1~X1, treat.2~X1, treat.3~X1))

multibin1<-CBPS(formulas.multi,method="over",type="MultiBin")

y<-cbind(treat.1,treat.2,treat.3) \%*\% c(2,2,2)+X1 \%*\% c(-2,8,7,6,2) + rnorm(n,sd=5)

contrast.mat<-rbind(
diag(3),
c(1,1,1),
c(1,1,-1),
c(1,-1,-1),
c(1,-1,1)
)

treat.mat<-cbind(treat.1,treat.2,treat.3)\%*\%t(contrast.mat)


Rsq.raw<-apply(treat.mat,2,F=function(x) summary(lm(x~X1))$r.sq)
Rsq.wt<-apply(treat.mat,2,F=function(x) summary(lm(x~X1,w=multibin1$w))$r.sq)

plot(Rsq.raw,1:7,type="n",xlim=range(cbind(Rsq.raw,Rsq.wt)),ylab="Contrasts",xlab="R-Squared")
abline(h=c(1:7),col=gray(.5))
points(Rsq.raw,1:7,pch="O")
points(Rsq.wt,1:7,pch="X")


#Unweighted model
lm1<-lm(y~treat.1+treat.2+treat.3+X1)
summary(lm1)


#CBMSM weighted model
lm2<-lm(y~treat.1+treat.2+treat.3+X1, weights=multibin1$weights)
summary(lm2)

}
}